from __future__ import annotations

import os
import re
from typing import Any, Dict, Optional


def check_parameters(
    mut: str,
    bam: str,
    fasta: str,
    sample_id: Optional[str],
    neg_offset_mate_search: int,
    pos_offset_mate_search: int,
    one_based: bool,
    flag_bam_list: Dict[str, Any],
    report_bam_info: bool,
    report_softclip: bool,
    report_5p_3p_bases_fragment: int,
    remove_softclip: bool,
    retain_fail_qc: bool,
    tmp_folder: Optional[str],
    output_path: Optional[str],
    verbose: bool,
    n_cores: int,
) -> None:
    """Validate all parameters before the analysis pipeline starts.

    Mirrors check_parameters() in R/check_parameters.R.  Raises ``ValueError``
    or ``FileNotFoundError`` on the first invalid parameter encountered.
    """
    _check_mut(mut)
    _check_bam(bam, verbose)
    _check_fasta(fasta, verbose)
    _check_sample(sample_id)
    _check_neg_offset(neg_offset_mate_search)
    _check_pos_offset(pos_offset_mate_search)
    _check_bool("one_based", one_based)
    _check_flag_bam_list(flag_bam_list)
    _check_bool("report_bam_info", report_bam_info)
    _check_bool("report_softclip", report_softclip)
    _check_nonneg_int("report_5p_3p_bases_fragment", report_5p_3p_bases_fragment)
    _check_bool("remove_softclip", remove_softclip)
    _check_bool("retain_fail_qc", retain_fail_qc)
    _check_tmp_folder(tmp_folder)
    _check_output_path(output_path)
    _check_bool("verbose", verbose)
    _check_n_cores(n_cores)


# ---------------------------------------------------------------------------
# Individual parameter checks
# ---------------------------------------------------------------------------

def _check_mut(mut: str) -> None:
    if isinstance(mut, os.PathLike):
        mut = os.fspath(mut)
    if not isinstance(mut, str) or not mut:
        raise ValueError("'mut' must be a non-empty string.")
    is_file_format = bool(re.search(r"\.(vcf|tsv)(\.gz)?$", mut))
    if is_file_format and not os.path.isfile(mut):
        raise FileNotFoundError(f"Mutation file does not exist: {mut}")


def _check_bam(bam: str, verbose: bool) -> None:
    if isinstance(bam, os.PathLike):
        bam = os.fspath(bam)
    if not isinstance(bam, str):
        raise ValueError("'bam' must be a string path.")
    if not bam.endswith(".bam"):
        raise ValueError(f"BAM file does not have a .bam extension: {bam}")
    if not os.path.isfile(bam):
        raise FileNotFoundError(f"BAM file does not exist: {bam}")

    bai = bam + ".bai"
    if not os.path.isfile(bai):
        if verbose:
            print("Creating BAM index…")
        try:
            import pysam  # type: ignore
            pysam.index(bam)
        except ImportError:
            # pysam not available (Windows) — delegate to samtools subprocess
            import subprocess
            result = subprocess.run(
                ["samtools", "index", bam],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"Could not create BAM index for '{bam}'. "
                    f"Please run 'samtools index {bam}' manually."
                )
        except Exception as exc:
            raise RuntimeError(
                f"Could not create BAM index for '{bam}'. "
                f"Please run 'samtools index {bam}' manually."
            ) from exc


def _check_fasta(fasta: str, verbose: bool) -> None:
    if isinstance(fasta, os.PathLike):
        fasta = os.fspath(fasta)
    if not isinstance(fasta, str):
        raise ValueError("'fasta' must be a string path.")
    if not re.search(r"\.fa(sta)?$", fasta):
        raise ValueError(
            f"FASTA file does not have a .fa or .fasta extension: {fasta}"
        )
    if not os.path.isfile(fasta):
        raise FileNotFoundError(f"FASTA file does not exist: {fasta}")
    # pyfaidx auto-creates .fai when the Fasta handle is opened.


def _check_sample(sample_id: Optional[str]) -> None:
    if sample_id is not None and sample_id == "":
        raise ValueError("'sample_id' cannot be an empty string. Use None instead.")


def _check_neg_offset(value: int) -> None:
    if not isinstance(value, int):
        raise TypeError("'neg_offset_mate_search' must be an integer.")
    if value > 0:
        raise ValueError("'neg_offset_mate_search' must be <= 0.")


def _check_pos_offset(value: int) -> None:
    if not isinstance(value, int):
        raise TypeError("'pos_offset_mate_search' must be an integer.")
    if value < 0:
        raise ValueError("'pos_offset_mate_search' must be >= 0.")


def _check_bool(name: str, value: Any) -> None:
    if not isinstance(value, bool):
        raise TypeError(f"'{name}' must be a boolean (True/False).")


def _check_nonneg_int(name: str, value: Any) -> None:
    if not isinstance(value, int) or value < 0:
        raise ValueError(f"'{name}' must be a non-negative integer.")


_VALID_FLAG_NAMES = {
    "isPaired", "isProperPair", "isUnmappedQuery", "hasUnmappedMate",
    "isMinusStrand", "isMateMinusStrand", "isFirstMateRead", "isSecondMateRead",
    "isSecondaryAlignment", "isSupplementaryAlignment",
    "isNotPassingQualityControls", "isDuplicate",
}


def _check_flag_bam_list(flag_bam_list: Any) -> None:
    if not isinstance(flag_bam_list, dict):
        raise TypeError("'flag_bam_list' must be a dict.")
    for key, val in flag_bam_list.items():
        if key not in _VALID_FLAG_NAMES:
            raise ValueError(
                f"Invalid flag name '{key}' in 'flag_bam_list'. "
                f"Valid names: {sorted(_VALID_FLAG_NAMES)}"
            )
        if val is not None and not isinstance(val, bool):
            raise TypeError(
                f"All values in 'flag_bam_list' must be bool or None, "
                f"got {type(val)} for key '{key}'."
            )


def _check_tmp_folder(tmp_folder: Optional[str]) -> None:
    if tmp_folder is None:
        return
    if not isinstance(tmp_folder, str):
        raise TypeError("'tmp_folder' must be a string or None.")
    os.makedirs(tmp_folder, exist_ok=True)


def _check_output_path(output_path: Optional[str]) -> None:
    if output_path is None or output_path == "":
        return
    if not isinstance(output_path, str):
        raise TypeError("'output_path' must be a string or None.")


def _check_n_cores(n_cores: int) -> None:
    if not isinstance(n_cores, int) or n_cores <= 0:
        raise ValueError("'n_cores' must be a positive integer.")
