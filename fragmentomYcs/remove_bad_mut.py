from __future__ import annotations

import re
import warnings
from typing import Any, Dict, List


def remove_bad_mut(df_mut: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Filter a list of variant dicts, keeping only rows that pass all validity
    checks.  Mirrors remove_bad_mut() in R/remove_bad_mut.R.

    Parameters
    ----------
    df_mut:
        List of dicts with keys ``CHROM``, ``POS``, ``REF``, ``ALT``.

    Returns
    -------
    list[dict]
        Filtered list containing only valid rows.

    Raises
    ------
    ValueError
        If no valid mutations remain after filtering.
    """
    valid: List[Dict[str, Any]] = []

    for i, row in enumerate(df_mut):
        chrom = row.get("CHROM")
        pos = row.get("POS")
        ref = row.get("REF")
        alt = row.get("ALT")

        if (
            _check_chr_input(chrom)
            and _check_pos_input(pos)
            and _check_ref_alt_input(ref, alt)
        ):
            valid.append(row)
        else:
            warnings.warn(
                f"Invalid row: {i + 1} - CHROM: {chrom} POS: {pos} REF: {ref} ALT: {alt}",
                stacklevel=2,
            )

    if not valid:
        raise ValueError("No valid mutations found after sanity check.")

    return valid


# ---------------------------------------------------------------------------
# Helper validators (mirrors R/remove_bad_mut.R)
# ---------------------------------------------------------------------------

def _check_chr_input(chr_: Any) -> bool:
    """Return True if *chr_* looks like a valid chromosome name."""
    if chr_ is None or chr_ == "":
        return False
    if isinstance(chr_, float) and chr_ != chr_:  # float NaN is not equal to itself
        return False
    return bool(
        re.match(
            r"^(chr([1-9]|1[0-9]|2[0-2]|X|Y)|([1-9]|1[0-9]|2[0-2]|X|Y))$",
            str(chr_),
        )
    )


def _check_pos_input(pos: Any) -> bool:
    """Return True if *pos* is a positive integer-valued number."""
    if pos is None or pos == "":
        return False
    if isinstance(pos, float) and pos != pos:  # NaN check
        return False
    try:
        p = float(pos)
        return p == int(p) and int(p) > 0
    except (TypeError, ValueError):
        return False


def _check_ref_alt_input(ref: Any, alt: Any) -> bool:
    """Return True if the ref/alt pair is a valid (possibly pre-VCF) allele pair.

    Mirrors check_ref_alt_input() in R: accepts ATCG characters as well as the
    placeholder characters ``.``, ``-``, ``_`` that are later normalised away.
    Rejects pairs where BOTH alleles are empty/missing, and rejects multi-allelic
    commas.
    """
    _specific = {"", ".", "-", "_", None, "NA", "nan"}

    ref_s = str(ref).strip() if ref is not None else ""
    alt_s = str(alt).strip() if alt is not None else ""

    # Treat Python NaN strings as missing
    if ref_s.lower() == "nan":
        ref_s = ""
    if alt_s.lower() == "nan":
        alt_s = ""

    # Both alleles are empty/placeholder → invalid
    if ref_s in _specific and alt_s in _specific:
        return False

    # Multi-allelic representations are not handled here
    if "," in ref_s or "," in alt_s:
        return False

    # Valid pattern: ATCG letters, optionally with one trailing placeholder
    valid_pattern = re.compile(r"^[ATCGatcg]*[._-]?$")
    return bool(valid_pattern.match(ref_s)) and bool(valid_pattern.match(alt_s))
