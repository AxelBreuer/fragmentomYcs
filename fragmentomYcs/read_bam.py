from __future__ import annotations

from typing import Any, Dict, List, Optional

try:
    import pysam as _bam_lib
except ImportError:
    try:
        import bamnostic as _bam_lib  # type: ignore[no-redef]
    except ImportError as exc:
        raise ImportError(
            "Neither pysam nor bamnostic is installed.\n"
            "  Linux/macOS: conda install -c conda-forge pysam\n"
            "  Windows:     pip install bamnostic"
        ) from exc


def read_bam(
    bam: str,
    chr_: str,
    pos: int,
    neg_offset_mate_search: int,
    pos_offset_mate_search: int,
    flag_bam_list: Dict[str, Any],
) -> Optional[List[Dict[str, Any]]]:
    """Extract and filter paired-end reads for a target locus from a BAM file.

    Mirrors the behaviour of read_bam() in R/read_bam.R:
    1. Defines an extended region around *pos* using the two offset parameters.
    2. Fetches all reads in that region that pass the FLAG filters in
       *flag_bam_list*.
    3. Identifies the subset of reads that directly cover *pos*.
    4. Collects the fragment names (QNAME) of those covering reads.
    5. Returns ALL reads from the initial fetch whose QNAME is in that set
       (i.e. mates are retrieved even if they do not cover *pos*).

    Parameters
    ----------
    bam:
        Path to an indexed BAM file.
    chr_:
        Chromosome of interest (must match the contig name in the BAM header).
    pos:
        1-based genomic position of interest.
    neg_offset_mate_search:
        Negative offset applied to *pos* to define the left boundary of the
        search window (pass a negative integer, e.g. -500).
    pos_offset_mate_search:
        Positive offset applied to *pos* to define the right boundary.
    flag_bam_list:
        Keyword arguments forwarded to pysam's view filters.  Recognised keys:
        ``is_paired``, ``is_proper_pair``, ``is_unmapped``,
        ``is_mate_unmapped``, ``is_reverse``, ``is_mate_reverse``,
        ``is_read1``, ``is_read2``, ``is_secondary``, ``is_qcfail``,
        ``is_duplicate``, ``is_supplementary``.
        Set a key to ``True`` to *require* the flag, ``False`` to *exclude* it.

    Returns
    -------
    list[dict] | None
        A list of dicts (one per read) with columns QNAME, FLAG, RNAME, POS,
        TLEN, MAPQ, CIGAR, RNEXT, PNEXT, SEQ, QUAL — or None if no reads
        cover the position of interest.
    """
    # --------------------------------------- Define region
    start_ext = max(1, pos + neg_offset_mate_search)
    end_ext = pos + pos_offset_mate_search

    # pysam uses 0-based half-open intervals
    start_0 = start_ext - 1
    end_0 = end_ext  # exclusive

    # --------------------------------------- Open BAM and fetch reads
    with _bam_lib.AlignmentFile(bam, "rb") as bam_file:
        reads = list(bam_file.fetch(chr_, start_0, end_0))

    if not reads:
        return None

    # --------------------------------------- Apply FLAG filters
    reads = _apply_flag_filters(reads, flag_bam_list)

    if not reads:
        return None

    # --------------------------------------- Convert to dicts
    df_sam_filtered = [_read_to_dict(r) for r in reads]

    # --------------------------------------- Find reads that cover *pos*
    # (1-based inclusive: read_start <= pos <= read_end)
    cover_position = [
        r["POS"] <= pos <= _read_end(r)
        for r in df_sam_filtered
    ]

    fragments_of_interest = set(
        r["QNAME"]
        for r, covers in zip(df_sam_filtered, cover_position)
        if covers
    )

    if not fragments_of_interest:
        return None

    # --------------------------------------- Select all reads belonging to
    # the fragments of interest (including mates)
    df_reads_of_covering_fragments = [
        r for r in df_sam_filtered if r["QNAME"] in fragments_of_interest
    ]

    return df_reads_of_covering_fragments


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Map Rsamtools/R flag names (used in flag_bam_list) → pysam attribute names
_RSAMTOOLS_TO_PYSAM: Dict[str, str] = {
    "isPaired":                    "is_paired",
    "isProperPair":                "is_proper_pair",
    "isUnmappedQuery":             "is_unmapped",
    "hasUnmappedMate":             "mate_is_unmapped",
    "isMinusStrand":               "is_reverse",
    "isMateMinusStrand":           "mate_is_reverse",
    "isFirstMateRead":             "is_read1",
    "isSecondMateRead":            "is_read2",
    "isSecondaryAlignment":        "is_secondary",
    "isSupplementaryAlignment":    "is_supplementary",
    "isNotPassingQualityControls": "is_qcfail",
    "isDuplicate":                 "is_duplicate",
}


def _apply_flag_filters(
    reads: List[Any],
    flag_bam_list: Dict[str, Any],
) -> List[Any]:
    """Filter reads using flag attributes.

    Each key in *flag_bam_list* maps to an AlignedSegment attribute name
    (e.g. ``is_paired``).  ``True`` = flag must be set; ``False`` = must not
    be set; ``None`` = skip this filter.
    """
    if not flag_bam_list:
        return reads

    filtered = []
    for read in reads:
        keep = True
        for r_name, required in flag_bam_list.items():
            if required is None:          # None means "don't filter on this flag"
                continue
            pysam_attr = _RSAMTOOLS_TO_PYSAM.get(r_name, r_name)
            if getattr(read, pysam_attr, None) != required:
                keep = False
                break
        if keep:
            filtered.append(read)

    return filtered


def _read_to_dict(read: Any) -> Dict[str, Any]:
    """Convert a pysam AlignedSegment to a plain dict matching the R column
    names used throughout fragmentomYcs (QNAME, FLAG, RNAME, POS, …)."""
    # pysam stores positions as 0-based; convert to 1-based to match R / SAM
    pos = (read.reference_start + 1) if read.reference_start is not None else None
    pnext = (read.next_reference_start + 1) if read.next_reference_start is not None else None

    rnext = read.next_reference_name
    if rnext is None:
        rnext = "*"
    elif rnext == read.reference_name:
        rnext = "="

    return {
        "QNAME": read.query_name,
        "FLAG": int(read.flag),
        "RNAME": read.reference_name,
        "POS": pos,
        "TLEN": getattr(read, "template_length", None) or getattr(read, "tlen", None),
        "MAPQ": read.mapping_quality,
        "CIGAR": read.cigarstring if read.cigarstring else "*",
        "RNEXT": rnext,
        "PNEXT": pnext,
        "SEQ": read.query_sequence,
        "QUAL": _qual_to_str(read.query_qualities),
    }


def _qual_to_str(qualities) -> str:
    """Convert pysam quality array to ASCII QUAL string (Phred+33)."""
    if qualities is None:
        return "*"
    return "".join(chr(q + 33) for q in qualities)


def _cigar_ref_width(cigar: str) -> int:
    """Return the number of reference bases consumed by a CIGAR string.

    Mirrors get_cigar_width() in R/read_bam.R.
    Operations that consume the reference: M D N = X
    """
    import re
    return sum(
        int(length)
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if op in ("M", "D", "N", "=", "X")
    )


def _read_end(r: Dict[str, Any]) -> int:
    """Return the last reference position (1-based inclusive) of a read."""
    if r["POS"] is None or r["CIGAR"] in (None, "*"):
        return r["POS"] or 0
    return r["POS"] + _cigar_ref_width(r["CIGAR"]) - 1
