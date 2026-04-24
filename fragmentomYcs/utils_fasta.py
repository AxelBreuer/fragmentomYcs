from __future__ import annotations

from typing import Any, Mapping, Optional


def get_seq_from_fasta(
    chr_: str,
    start: int,
    end: int,
    fasta_seq: Optional[Mapping[str, Any]] = None,
) -> str:
    if fasta_seq is None:
        raise ValueError("fasta_seq is required in this Python implementation")

    f_chr = str(fasta_seq["chr"])
    f_start = int(fasta_seq["start"])
    f_end = int(fasta_seq["end"])
    f_seq = str(fasta_seq["seq"])

    if f_chr != chr_:
        raise ValueError(
            f"Requested chromosome '{chr_}' does not match available chromosome '{f_chr}'"
        )
    if start < f_start or end > f_end:
        raise ValueError(
            f"Requested sequence range {start}:{end} does not fit available range {f_start}:{f_end}"
        )

    start_idx = start - f_start
    end_idx = end - f_start + 1
    return f_seq[start_idx:end_idx]
