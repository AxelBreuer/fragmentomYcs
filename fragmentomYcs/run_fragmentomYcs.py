"""Main entry point for the fragmentomYcs Python pipeline.

Mirrors the R ``run_fragmentomYcs`` function.  The normalisation backend
(pybcftools or system bcftools binary) is configured via ``fragmentomYcs.cfg``.

Parallel processing uses :class:`concurrent.futures.ProcessPoolExecutor`.

.. note::
    On Windows, multiprocessing requires caller scripts to use an
    ``if __name__ == "__main__":`` guard.  Single-core mode (``n_cores=1``)
    needs no guard.
"""

from __future__ import annotations

import os
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

import pandas as pd
import pyfaidx

from fragmentomYcs.check_parameters import check_parameters
from fragmentomYcs.extract_fragment_features import extract_fragment_features
from fragmentomYcs.normalize_mut import normalize_mut
from fragmentomYcs.normalize_to_vcf_rep import get_seq_from_fasta_file
from fragmentomYcs.read_bam import read_bam
from fragmentomYcs.read_mut import read_mut
from fragmentomYcs.remove_bad_mut import remove_bad_mut


# Number of fragment sub-batches dispatched to workers per variant
_N_CHUNKS = 50


def _chunk(lst: list, size: int) -> list:
    """Split *lst* into consecutive sub-lists of at most *size* elements."""
    return [lst[k : k + size] for k in range(0, len(lst), size)]


# ---------------------------------------------------------------------------
# Module-level worker (must be at module level to be picklable on Windows)
# ---------------------------------------------------------------------------

def _process_chunk(args: tuple) -> List[Dict[str, Any]]:
    """Process one chunk of fragments.  Called by worker processes."""
    (
        df_sam_chunk,
        fragment_chunk,
        sample_id,
        chr_norm,
        pos_norm,
        ref_norm,
        alt_norm,
        report_bam_info,
        report_softclip,
        report_5p_3p_bases_fragment,
        remove_softclip_flag,
        fasta_seq,
        input_mutation_info,
    ) = args

    results: List[Dict[str, Any]] = []
    for fragment_name in fragment_chunk:
        result = extract_fragment_features(
            df_sam=df_sam_chunk,
            fragment_name=fragment_name,
            sample_id=sample_id,
            chr_=chr_norm,
            pos=pos_norm,
            ref=ref_norm,
            alt=alt_norm,
            report_bam_info=report_bam_info,
            report_softclip=report_softclip,
            report_5p_3p_bases_fragment=report_5p_3p_bases_fragment,
            remove_softclip=remove_softclip_flag,
            fasta_seq=fasta_seq,
            input_mutation_info=input_mutation_info,
        )
        results.append(result)
    return results


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_fragmentomYcs(
    mut: str,
    bam: str,
    fasta: str,
    sample_id: Optional[str] = None,
    neg_offset_mate_search: int = -600,
    pos_offset_mate_search: int = 600,
    one_based: bool = True,
    flag_bam_list: Optional[Dict[str, Any]] = None,
    report_bam_info: bool = False,
    report_softclip: bool = False,
    report_5p_3p_bases_fragment: int = 5,
    remove_softclip: bool = False,
    retain_fail_qc: bool = False,
    apply_bcftools_norm: bool = False,
    tmp_folder: Optional[str] = None,
    output_path: Optional[str] = None,
    verbose: bool = False,
    n_cores: int = 1,
) -> Optional[pd.DataFrame]:
    """End-to-end pipeline for fragment-level allelic analysis.

    Parameters
    ----------
    mut:
        Path to a ``.vcf`` / ``.tsv`` mutation file, or a
        ``"chr:pos:ref:alt"`` string.
    bam:
        Path to an indexed BAM file.
    fasta:
        Path to the reference FASTA file used to generate the BAM.
    sample_id:
        Optional sample identifier added to every output row.
    neg_offset_mate_search:
        Upstream search window offset (negative integer, default -600).
    pos_offset_mate_search:
        Downstream search window offset (positive integer, default 600).
    one_based:
        True if positions in *mut* are 1-based (default True).
    flag_bam_list:
        SAM-flag filter dict.  Keys are Rsamtools-style flag names (e.g.
        ``"isPaired"``); values are ``True`` (flag must be set), ``False``
        (flag must be unset), or ``None`` (ignored).
    report_bam_info:
        Include raw BAM fields in the output.
    report_softclip:
        Include soft-clip counts at fragment extremities.
    report_5p_3p_bases_fragment:
        Number of fragment-end bases to report (default 5).
    remove_softclip:
        Trim soft-clipped bases before analysis.
    retain_fail_qc:
        Keep fragments that failed QC in the output (default False).
    apply_bcftools_norm:
        If True, left-align each variant. The backend (pybcftools or bcftools
        binary) is selected via ``fragmentomYcs.cfg``.
    tmp_folder:
        Temporary directory path (defaults to ``tempfile.gettempdir()``).
    output_path:
        If provided, results are written to this tab-separated file and
        the function returns ``None``.
    verbose:
        Print progress messages.
    n_cores:
        Number of parallel worker processes (default 1).

    Returns
    -------
    pandas.DataFrame or None
        DataFrame with one row per fragment, or ``None`` if *output_path*
        was specified.

    Raises
    ------
    ValueError
        If any input parameter is invalid or no reads were found.
    """
    # ------------------------------------------------------------------
    # 0. Default mutable argument
    # ------------------------------------------------------------------
    if flag_bam_list is None:
        flag_bam_list = {
            "isPaired": True,
            "isProperPair": None,
            "isUnmappedQuery": False,
            "hasUnmappedMate": False,
            "isMinusStrand": None,
            "isMateMinusStrand": None,
            "isFirstMateRead": None,
            "isSecondMateRead": None,
            "isSecondaryAlignment": False,
            "isSupplementaryAlignment": False,
            "isNotPassingQualityControls": None,
            "isDuplicate": None,
        }

    if tmp_folder is None:
        tmp_folder = tempfile.gettempdir()

    # Normalise path-like objects to plain strings for all downstream code
    import os as _os
    mut   = _os.fspath(mut)
    bam   = _os.fspath(bam)
    fasta = _os.fspath(fasta)

    # ------------------------------------------------------------------
    # 1. Validate parameters
    # ------------------------------------------------------------------
    check_parameters(
        mut=mut,
        bam=bam,
        fasta=fasta,
        sample_id=sample_id,
        neg_offset_mate_search=neg_offset_mate_search,
        pos_offset_mate_search=pos_offset_mate_search,
        one_based=one_based,
        flag_bam_list=flag_bam_list,
        report_bam_info=report_bam_info,
        report_softclip=report_softclip,
        report_5p_3p_bases_fragment=report_5p_3p_bases_fragment,
        remove_softclip=remove_softclip,
        retain_fail_qc=retain_fail_qc,
        tmp_folder=tmp_folder,
        output_path=output_path,
        verbose=verbose,
        n_cores=n_cores,
    )

    # ------------------------------------------------------------------
    # 2. Load mutations and open FASTA
    # ------------------------------------------------------------------
    df_mut_raw = read_mut(mut)
    df_mut_raw = remove_bad_mut(df_mut_raw)

    # pyfaidx automatically creates the .fai index if missing
    fa = pyfaidx.Fasta(fasta)

    df_mut_norm = normalize_mut(
        df_mut=df_mut_raw,
        fasta=fa,
        fasta_path=fasta,
        one_based=one_based,
        apply_bcftools_norm=apply_bcftools_norm,
        verbose=verbose,
    )

    if not df_mut_norm:
        raise ValueError("No mutations remained after normalisation.")

    # ------------------------------------------------------------------
    # 3. Per-mutation analysis
    # ------------------------------------------------------------------
    final_n_cores = min(n_cores, os.cpu_count() or 1)
    results_list: List[Optional[pd.DataFrame]] = [None] * len(df_mut_norm)

    for i, mut_row in enumerate(df_mut_norm):
        chr_norm: str = mut_row["chr"]
        pos_norm: int = int(mut_row["pos"])
        ref_norm: str = mut_row["ref"]
        alt_norm: str = mut_row["alt"]
        input_mutation_info: str = mut_row["input_mutation_info"]

        if verbose:
            print(
                f"[{i + 1}/{len(df_mut_norm)}] Processing "
                f"{chr_norm}:{pos_norm} {ref_norm}>{alt_norm} "
                f"(input: {input_mutation_info})"
            )

        # ---- 3a. Extract reads from BAM
        df_sam = read_bam(
            bam=bam,
            chr_=chr_norm,
            pos=pos_norm,
            neg_offset_mate_search=neg_offset_mate_search,
            pos_offset_mate_search=pos_offset_mate_search,
            flag_bam_list=flag_bam_list,
        )

        if df_sam is None:
            warnings.warn(
                f"No reads cover {chr_norm}:{pos_norm}:{ref_norm}>{alt_norm}. "
                "Skipping.",
                stacklevel=2,
            )
            continue

        # ---- 3b. Pre-load FASTA window for this variant
        #          (avoids per-fragment FASTA I/O inside workers)
        if len(ref_norm) == len(alt_norm):
            motif_len = len(alt_norm)
        else:
            motif_len = max(len(ref_norm) - 1, len(alt_norm) - 1)

        seq_lengths = [len(r["SEQ"]) for r in df_sam if r.get("SEQ")]
        max_read_seq_len = max(seq_lengths) if seq_lengths else 150

        same_chr = [r for r in df_sam if r.get("RNAME") == chr_norm]
        if not same_chr:
            warnings.warn(
                f"No reads on chromosome {chr_norm} for position {pos_norm}.",
                stacklevel=2,
            )
            continue

        pos_vals = [r["POS"] for r in same_chr if r.get("POS") is not None]
        min_read_pos = min(pos_vals)

        min_fetch_pos = min_read_pos - 1  # extra base for look-ahead
        max_fetch_pos = max(pos_vals) + max_read_seq_len + motif_len + 1

        fetch_seq = get_seq_from_fasta_file(
            fa, chr_norm, min_fetch_pos, max_fetch_pos
        )
        fasta_seq: Dict[str, Any] = {
            "chr": chr_norm,
            "start": min_fetch_pos,
            "end": max_fetch_pos,
            "seq": fetch_seq,
        }

        # ---- 3c. Parallel fragment processing (chunked)
        fragment_names = list({r["QNAME"] for r in df_sam})

        chunk_size = max(1, -(-len(fragment_names) // _N_CHUNKS))  # ceil division
        fragment_chunks = _chunk(fragment_names, chunk_size)

        # Pre-split reads by chunk to reduce data sent to each worker
        df_sam_chunks = [
            [r for r in df_sam if r["QNAME"] in set(ch)]
            for ch in fragment_chunks
        ]

        chunk_args = [
            (
                df_sam_chunks[idx],
                fragment_chunks[idx],
                sample_id,
                chr_norm,
                pos_norm,
                ref_norm,
                alt_norm,
                report_bam_info,
                report_softclip,
                report_5p_3p_bases_fragment,
                remove_softclip,
                fasta_seq,
                input_mutation_info,
            )
            for idx in range(len(fragment_chunks))
        ]

        all_fragment_results: List[Dict[str, Any]] = []

        if final_n_cores > 1:
            with ProcessPoolExecutor(max_workers=final_n_cores) as executor:
                futures = {
                    executor.submit(_process_chunk, arg): arg for arg in chunk_args
                }
                for future in as_completed(futures):
                    all_fragment_results.extend(future.result())
        else:
            for arg in chunk_args:
                all_fragment_results.extend(_process_chunk(arg))

        if not all_fragment_results:
            continue

        df_frags = pd.DataFrame(all_fragment_results)

        # ---- 3d. Compute VAF
        status_simple = df_frags.get("Fragment_Status_Simple", pd.Series(dtype=str))
        if status_simple.isna().all():
            df_frags["VAF"] = float("nan")
        else:
            denom = status_simple.isin(["MUT", "WT", "OTH"]).sum()
            num = (status_simple == "MUT").sum()
            df_frags["VAF"] = 0.0 if denom == 0 else 100.0 * num / denom

        results_list[i] = df_frags

    # ------------------------------------------------------------------
    # 4. Combine all per-mutation results
    # ------------------------------------------------------------------
    valid_dfs = [df for df in results_list if df is not None]
    if not valid_dfs:
        raise ValueError("The final fragmentomYcs dataframe is empty.")

    df_final = pd.concat(valid_dfs, ignore_index=True)

    # ------------------------------------------------------------------
    # 5. Remove QC-failed fragments (unless requested otherwise)
    # ------------------------------------------------------------------
    if not retain_fail_qc:
        mask_fail = df_final["Fragment_QC"] != "OK"
        n_fail = int(mask_fail.sum())
        df_final = df_final.loc[~mask_fail].reset_index(drop=True)
        if verbose and n_fail:
            print(
                f"Removed {n_fail} fragment(s) that failed quality checks. "
                "Set retain_fail_qc=True to retain them."
            )

    # ------------------------------------------------------------------
    # 6. Write output file or return DataFrame
    # ------------------------------------------------------------------
    if output_path:
        parent = os.path.dirname(output_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        if os.path.isfile(output_path) and verbose:
            print(f"File '{output_path}' already exists and will be overwritten.")
        if verbose:
            print(f"Writing results to: {output_path}")
        df_final.to_csv(output_path, sep="\t", index=False)
        return None

    return df_final
