from __future__ import annotations

from typing import Any, Mapping

from fragmentomYcs.utils import parse_cigar


def get_index_aligning_with_pos(pos: int, read_stats: Mapping[str, Any]) -> float:
    read_pos = int(read_stats["POS"])
    cigar_ops = parse_cigar(str(read_stats["CIGAR"]))

    operations = {
        "M": (True, True),
        "I": (True, False),
        "D": (False, True),
        "N": (False, True),
        "S": (True, False),
        "H": (False, False),
        "P": (False, False),
        "=": (True, True),
        "X": (True, True),
    }

    ref_pos = read_pos
    read_idx = 1

    for op_idx, (op_len, op_type) in enumerate(cigar_ops):
        consumes_seq, consumes_ref = operations[op_type]

        ref_move = op_len if consumes_ref else 0
        seq_move = op_len if consumes_seq else 0

        if consumes_ref and pos >= ref_pos and pos < (ref_pos + ref_move):
            if op_type in ("M", "=", "X"):
                return float(read_idx + (pos - ref_pos))
            if op_type in ("D", "N"):
                return -2.0

        if consumes_seq and pos == ref_pos and op_idx == len(cigar_ops) - 1:
            return -0.5

        ref_pos += ref_move
        read_idx += seq_move

    return -1.0
