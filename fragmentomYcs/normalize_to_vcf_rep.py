from __future__ import annotations

import re
import warnings
from typing import Any, Dict, Optional

import pyfaidx


def normalize_to_vcf_rep(
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    fasta: pyfaidx.Fasta,
    one_based: bool,
) -> Optional[Dict[str, Any]]:
    """Convert a user-supplied variant to canonical VCF representation.

    Strips allele placeholders, adjusts *pos* to 1-based, harmonises the
    chromosome name to match the FASTA header, adds anchor bases for pure
    indels, and validates *ref* against the FASTA.
    Mirrors normalize_to_vcf_rep() in R/normalize_to_vcf_rep.R.

    Parameters
    ----------
    chr_, pos, ref, alt:
        Variant as read from the input file.
    fasta:
        Open ``pyfaidx.Fasta`` handle for the reference genome.
    one_based:
        True if *pos* is already 1-based (standard for VCF/TSV inputs).

    Returns
    -------
    dict with keys ``chr``, ``pos``, ``ref``, ``alt``, or ``None`` if *ref*
    does not match the FASTA sequence.
    """
    # ---- 1. Strip placeholder characters
    ref_norm, alt_norm = _normalize_na_representation(ref, alt)

    # ---- 2. Convert 0-based → 1-based if needed
    pos_norm: int = pos if one_based else pos + 1

    # ---- 3. Harmonise chromosome name
    chr_norm = harmonize_chr_to_fasta(str(chr_), fasta)

    # ---- 4. Add anchor base for pure indels (VCF convention)
    if len(ref_norm) != len(alt_norm):
        if len(ref_norm) > len(alt_norm):
            # Deletion: alt allele was empty → prepend anchor at pos_norm - 1
            if len(alt_norm) == 0:
                anchor = _fetch_base(fasta, chr_norm, pos_norm - 1)
                alt_norm = anchor
                ref_norm = anchor + ref_norm
                pos_norm -= 1
        else:
            # Insertion: ref allele was empty → prepend anchor at pos_norm
            if len(ref_norm) == 0:
                anchor = _fetch_base(fasta, chr_norm, pos_norm)
                ref_norm = anchor
                alt_norm = anchor + alt_norm

    # ---- 5. Validate ref against FASTA
    if not _check_if_ref_matches_fasta(chr_norm, pos_norm, ref_norm, fasta):
        warnings.warn(
            f"Mismatch between ref allele and FASTA for ({chr_} {pos} {ref}). "
            "Skipping the mutation.",
            stacklevel=2,
        )
        return None

    return {"chr": chr_norm, "pos": pos_norm, "ref": ref_norm, "alt": alt_norm}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _normalize_na_representation(ref: str, alt: str):
    """Strip placeholder characters (mirrors normalize_na_representation in R)."""
    ref = re.sub(r"[-._]|NA", "", str(ref))
    alt = re.sub(r"[-._]|NA", "", str(alt))
    return ref, alt


def harmonize_chr_to_fasta(chr_: str, fasta: pyfaidx.Fasta) -> str:
    """Adjust *chr_* prefix to match the style used in *fasta* headers."""
    fasta_chrs = list(fasta.keys())
    has_chr_prefix = any(k.startswith("chr") for k in fasta_chrs)

    if has_chr_prefix and not chr_.startswith("chr"):
        return "chr" + chr_
    if not has_chr_prefix and chr_.startswith("chr"):
        return chr_[3:]
    return chr_


def _fetch_base(fasta: pyfaidx.Fasta, chr_: str, pos: int) -> str:
    """Fetch a single nucleotide at 1-based *pos* from *fasta*."""
    # pyfaidx uses 0-based half-open: [pos-1, pos)
    return str(fasta[chr_][pos - 1 : pos])


def get_seq_from_fasta_file(
    fasta: pyfaidx.Fasta, chr_: str, start: int, end: int
) -> str:
    """Fetch the reference sequence for the 1-based inclusive range [start, end].

    This is used during the normalisation phase where we need direct FASTA
    access (as opposed to the pre-loaded window used inside workers).
    """
    return str(fasta[chr_][start - 1 : end])


def _check_if_ref_matches_fasta(
    chr_: str, pos: int, ref: str, fasta: pyfaidx.Fasta
) -> bool:
    """Return True if *ref* matches the FASTA sequence at [pos, pos+len(ref))."""
    if chr_ not in fasta:
        warnings.warn(f"Chromosome '{chr_}' not found in FASTA.", stacklevel=3)
        return False
    fasta_seq = get_seq_from_fasta_file(fasta, chr_, pos, pos + len(ref) - 1)
    return fasta_seq.upper() == ref.upper()
