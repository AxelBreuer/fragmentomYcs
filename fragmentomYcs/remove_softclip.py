from __future__ import annotations

import re
from typing import Any, Dict, Mapping


def remove_softclip(read_stats: Mapping[str, Any]) -> Dict[str, Any]:
    cigar = str(read_stats["CIGAR"])
    seq = str(read_stats["SEQ"])
    qual = str(read_stats["QUAL"])

    softclip_5p = re.match(r"^(?:\d+H)?(\d+)S", cigar)
    n_softclip_5p = int(softclip_5p.group(1)) if softclip_5p else 0

    softclip_3p = re.search(r"(\d+)S(?:\d+H)?$", cigar)
    n_softclip_3p = int(softclip_3p.group(1)) if softclip_3p else 0

    seq_len = len(seq)
    start = n_softclip_5p
    end = seq_len - n_softclip_3p

    new_seq = seq[start:end] if end >= start else ""
    new_qual = qual[start:end] if end >= start else ""

    new_cigar = re.sub(r"^(\d+H)?\d+S", lambda m: m.group(1) or "", cigar)
    new_cigar = re.sub(r"\d+S(\d+H)?$", lambda m: m.group(1) or "", new_cigar)

    return {
        "SEQ": new_seq,
        "QUAL": new_qual,
        "CIGAR": new_cigar,
        "read_length": len(new_seq),
    }
