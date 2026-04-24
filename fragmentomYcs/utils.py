from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Mapping, Sequence, Tuple


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    matches = re.findall(r"([0-9]+)([MIDNSHP=X])", cigar)
    return [(int(length), op_type) for length, op_type in matches]


def get_pos_indels_from_read(read_stats: Mapping[str, Any]) -> Dict[str, List[float]]:
    list_pos_del: List[float] = []
    list_pos_ins: List[float] = []

    cigar_ops = parse_cigar(str(read_stats["CIGAR"]))
    current_ref_pos = int(read_stats["POS"])

    for op_len, op_type in cigar_ops:
        if op_type in ("M", "=", "X"):
            current_ref_pos += op_len
        elif op_type == "D":
            list_pos_del.extend(float(current_ref_pos + i) for i in range(op_len))
            current_ref_pos += op_len
        elif op_type == "I":
            pos_before_insertion = current_ref_pos - 1
            scale = 10 ** (math.ceil(math.log10(op_len + 1)))
            list_pos_ins.extend(pos_before_insertion + ((i + 1) / scale) for i in range(op_len))
        elif op_type == "N":
            current_ref_pos += op_len

    return {"deletions": list_pos_del, "insertions": list_pos_ins}


def sam_flag_bits(flag: int) -> Dict[str, bool]:
    return {
        "isMinusStrand": bool(flag & 0x10),
        "isMateMinusStrand": bool(flag & 0x20),
        "isFirstMateRead": bool(flag & 0x40),
        "isSecondMateRead": bool(flag & 0x80),
    }


def return_fail_qc_fragment(
    qc_message: str,
    sample_id: str | None,
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    input_mutation_info: str,
    fragment_name: str,
) -> Dict[str, Any]:
    return {
        "Sample_Id": sample_id,
        "Chromosome": chr_,
        "Position": int(pos),
        "Ref": ref,
        "Alt": alt,
        "Input_Mutation": input_mutation_info,
        "Fragment_Id": fragment_name,
        "Fragment_QC": qc_message,
        "Fragment_Status_Simple": None,
        "Fragment_Status_Detail": None,
        "Fragment_Size": None,
        "Read_5p_Status": None,
        "Read_3p_Status": None,
        "BASE_5p": None,
        "BASE_3p": None,
        "BASQ_5p": None,
        "BASQ_3p": None,
        "Position_5p": None,
        "Position_3p": None,
        "VAF": None,
    }


def get_number_of_common_first_char(str_a: str, str_b: str) -> int:
    i = 0
    while i < len(str_a) and i < len(str_b) and str_a[i] == str_b[i]:
        i += 1
    return i


def is_na(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def as_records(df_sam: Any) -> List[Dict[str, Any]]:
    if hasattr(df_sam, "to_dict"):
        return [dict(x) for x in df_sam.to_dict(orient="records")]

    if isinstance(df_sam, Sequence):
        return [dict(x) for x in df_sam]

    raise TypeError("df_sam must be a pandas.DataFrame or a sequence of mappings")
