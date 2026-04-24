from __future__ import annotations

import re
from typing import Any, Mapping

from fragmentomYcs.utils import get_pos_indels_from_read


def get_fragment_size(read_stats_5p: Mapping[str, Any], read_stats_3p: Mapping[str, Any]) -> int:
    # --- Extract necessary metrics ---
    cigar_5p = str(read_stats_5p["CIGAR"])
    cigar_3p = str(read_stats_3p["CIGAR"])

    # bases matched read
    bases_match_5p = sum(int(x) for x in re.findall(r"(\d+)(?=M)", cigar_5p))

    # bases deleted read
    bases_del_5p = sum(int(x) for x in re.findall(r"(\d+)(?=D)", cigar_5p))

    # bases soft-clipped left
    m_soft_left_3p = re.search(r"^(\d+)(?=S)", cigar_3p)
    bases_softcl_left_3p = int(m_soft_left_3p.group(1)) if m_soft_left_3p else 0

    # bases soft-clipped right
    m_soft_right_5p = re.search(r"(\d+)(?=S$)", cigar_5p)
    bases_softcl_right_5p = int(m_soft_right_5p.group(1)) if m_soft_right_5p else 0

    # --- Check presence of insertions at 3p of the 5p read and at 5p of the 3p read ---
    # 5p of read 3p
    m_ins_left_3p = re.search(r"^(\d+)(?=I)", cigar_3p)
    bases_ins_left_3p = int(m_ins_left_3p.group(1)) if m_ins_left_3p else 0

    # 3p of read 5p
    m_ins_rigth_5p = re.search(r"(\d+)(?=I$)", cigar_5p)
    bases_ins_rigth_5p = int(m_ins_rigth_5p.group(1)) if m_ins_rigth_5p else 0

    # --- Define overlapping windows ---
    start_overlap = int(read_stats_3p["POS"]) - bases_softcl_left_3p
    end_overlap = int(read_stats_5p["POS"]) + bases_match_5p + bases_del_5p + bases_softcl_right_5p - 1

    # Define the real range of the overlapping with inserted bases
    extended_start_overlap = start_overlap - bases_ins_left_3p
    extended_end_overlap = end_overlap + bases_ins_rigth_5p

    # --- Calculate inner size of the fragment ---
    inner_distance = extended_start_overlap - extended_end_overlap - 1

    # Taking indels into account in the overlapping windows to not count them
    # twice. Get indels for each read
    indels_5p = get_pos_indels_from_read(read_stats_5p)
    indels_3p = get_pos_indels_from_read(read_stats_3p)

    # Merge the deletion lists and keep unique positions in the overlapping
    # windows
    unique_fragment_deletions = set(indels_5p["deletions"]).intersection(indels_3p["deletions"])
    overlapping_window_deletion = [
        x for x in unique_fragment_deletions if start_overlap <= x <= end_overlap
    ]

    # Merge the insertion lists and keep unique positions in the overlapping
    # windows
    unique_fragment_insertions = set(indels_5p["insertions"]).intersection(indels_3p["insertions"])
    overlapping_window_insertion = [
        x for x in unique_fragment_insertions if start_overlap <= x <= end_overlap
    ]

    fragment_size = (
        int(read_stats_5p["read_length"])
        + int(read_stats_3p["read_length"])
        + inner_distance
        + len(overlapping_window_deletion)
        - len(overlapping_window_insertion)
    )
    return int(fragment_size)
