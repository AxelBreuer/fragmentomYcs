from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional

import pyfaidx

from fragmentomYcs.apply_bcftools_norm import apply_bcftools_norm as _apply_bcftools_norm
from fragmentomYcs.normalize_to_vcf_rep import normalize_to_vcf_rep


def normalize_mut(
    df_mut: List[Dict[str, Any]],
    fasta: pyfaidx.Fasta,
    fasta_path: str,
    one_based: bool,
    apply_bcftools_norm: bool,
    verbose: bool,
) -> List[Dict[str, Any]]:
    """Normalise a list of variant dicts to canonical VCF representation.

    For each variant:
    1. ``normalize_to_vcf_rep`` converts user-supplied notation to VCF format
       (anchor-base padding, coordinate adjustment, chr harmonisation).
    2. Optionally, ``apply_bcftools_norm`` left-aligns indels.  The backend
       (pybcftools or system bcftools binary) is selected via
       ``fragmentomYcs.cfg`` (see ``apply_bcftools_norm`` module).

    Mirrors normalize_mut() in R/normalize_mut.R.

    Parameters
    ----------
    df_mut:
        List of dicts with keys ``CHROM``, ``POS``, ``REF``, ``ALT``.
    fasta:
        Open ``pyfaidx.Fasta`` handle for the reference genome.
    fasta_path:
        Filesystem path to the FASTA file.  Forwarded to the ``bcftools``
        binary backend when it is selected.
    one_based:
        True if *POS* values are already 1-based.
    apply_bcftools_norm:
        If True, further left-align each variant after VCF normalisation.
    verbose:
        Print per-variant progress messages.

    Returns
    -------
    list[dict]
        Each dict has keys ``chr``, ``pos``, ``ref``, ``alt``,
        ``input_mutation_info``.  Rows that fail validation are silently
        dropped (a warning is emitted for each).
    """
    results: List[Dict[str, Any]] = []

    for row in df_mut:
        chr_ = str(row["CHROM"])
        pos = int(row["POS"])
        ref = str(row["REF"])
        alt = str(row["ALT"])

        input_mutation_info = f"{chr_}:{pos}:{ref}-{alt}"

        if verbose:
            print(f"  Normalising {input_mutation_info}")

        # Step 1: VCF normalisation (anchor base, coordinate shift, chr name)
        vcf_norm: Optional[Dict[str, Any]] = normalize_to_vcf_rep(
            chr_=chr_,
            pos=pos,
            ref=ref,
            alt=alt,
            fasta=fasta,
            one_based=one_based,
        )

        if vcf_norm is None:
            warnings.warn(
                f"VCF normalisation failed for {input_mutation_info}. Skipping.",
                stacklevel=2,
            )
            continue

        # Step 2 (optional): left-align indels
        if apply_bcftools_norm:
            norm: Optional[Dict[str, Any]] = _apply_bcftools_norm(
                chr_=vcf_norm["chr"],
                pos=vcf_norm["pos"],
                ref=vcf_norm["ref"],
                alt=vcf_norm["alt"],
                fasta=fasta,
                fasta_path=fasta_path,
                verbose=verbose,
            )
            if norm is None:
                warnings.warn(
                    f"Left-alignment failed for {input_mutation_info}. Skipping.",
                    stacklevel=2,
                )
                continue
            norm["input_mutation_info"] = input_mutation_info
            results.append(norm)
        else:
            vcf_norm["input_mutation_info"] = input_mutation_info
            results.append(vcf_norm)

    return results
