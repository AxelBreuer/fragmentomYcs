from __future__ import annotations

from typing import Any, Mapping, Tuple

from fragmentomYcs.utils import parse_cigar


def search_for_indel_in_cigar(
    pos: int,
    ref: str,
    alt: str,
    read_stats: Mapping[str, Any],
    type_: str,
) -> Tuple[bool, bool]:
    if type_ not in ("I", "D"):
        raise ValueError("type_ must be 'I' or 'D'")

    ops = parse_cigar(str(read_stats["CIGAR"]))
    ref_pos = int(read_stats["POS"])
    read_idx = 0
    seq = str(read_stats["SEQ"])

    for op_len, op_type in ops:
        if op_type == type_:
            if op_type == "I" and ref_pos - 1 == pos:
                if op_len == len(alt) - 1:
                    left_base_idx = max(read_idx - 1, 0)
                    candidate = seq[left_base_idx : left_base_idx + op_len + 1]
                    if candidate == alt:
                        return True, False
                return False, True

            if op_type == "D" and ref_pos - 1 == pos:
                if op_len == len(ref) - 1:
                    return True, False
                return False, True

        if op_type in ("M", "D", "N", "=", "X"):
            ref_pos += op_len
        if op_type in ("M", "I", "S", "=", "X"):
            read_idx += op_len

    return False, False
