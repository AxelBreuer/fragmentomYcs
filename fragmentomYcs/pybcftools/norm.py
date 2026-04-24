"""Pure-Python left-alignment of VCF variants (Tan et al. 2015).

Implements the same normalisation performed by::

    bcftools norm -m +both -d exact --check REF,ALT -f <fasta>

without requiring the ``bcftools`` binary.  Works on Linux, macOS, and
Windows.

Reference
---------
Tan A, Abecasis GR, Kang HM.  Unified representation of genetic variants.
Bioinformatics. 2015;31(13):2202-4.  https://doi.org/10.1093/bioinformatics/btv112
"""

from __future__ import annotations

from typing import Tuple

import pyfaidx


def left_align_and_trim(
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    fasta: pyfaidx.Fasta,
) -> Tuple[int, str, str]:
    """Left-align and minimally represent a biallelic VCF variant.

    Parameters
    ----------
    chr_ : str
        Chromosome name (must exist in *fasta*).
    pos : int
        1-based genomic position (VCF convention).
    ref : str
        Reference allele in VCF format (anchor base included for indels).
    alt : str
        Alternate allele in VCF format (anchor base included for indels).
    fasta : pyfaidx.Fasta
        Open reference FASTA handle.

    Returns
    -------
    tuple[int, str, str]
        ``(pos, ref, alt)`` after left-alignment and trimming.

    Algorithm (Tan et al. 2015)
    ---------------------------
    Deletion:  shift left while the base immediately after the span equals
               the anchor base.  Each step prepends a new anchor and drops
               the rightmost base of *ref*.
    Insertion: shift left while the last inserted base equals the anchor.
               Each step prepends a new anchor and drops the rightmost base
               of *alt*.
    SNV/MNV:   trim identical trailing bases, then identical leading bases.
    """
    is_deletion = len(ref) > len(alt)
    is_insertion = len(alt) > len(ref)

    # ---- SNV / MNV: no left-alignment, just trim common flanks
    if not (is_deletion or is_insertion):
        while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
        return pos, ref, alt

    # ---- Deletion: shift left while base-after-deletion == anchor
    if is_deletion:
        while pos > 1:
            # 1-based position of the base immediately after the deletion:
            #   pos + len(ref) - 1
            # As a 0-based pyfaidx index:
            #   pos + len(ref) - 2
            base_after = _fetch(fasta, chr_, pos + len(ref) - 2)
            if ref[0] != base_after:
                break
            new_anchor = _fetch(fasta, chr_, pos - 2)  # genome[pos-1] in 0-based
            ref = new_anchor + ref[:-1]
            alt = new_anchor
            pos -= 1

    # ---- Insertion: shift left while last inserted base == anchor
    else:
        while pos > 1:
            if ref[0] != alt[-1]:
                break
            new_anchor = _fetch(fasta, chr_, pos - 2)  # genome[pos-1] in 0-based
            alt = new_anchor + alt[:-1]
            ref = new_anchor
            pos -= 1

    return pos, ref, alt


# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

def _fetch(fasta: pyfaidx.Fasta, chr_: str, pos0: int) -> str:
    """Return the single nucleotide at 0-based index *pos0*."""
    return str(fasta[chr_][pos0 : pos0 + 1])
