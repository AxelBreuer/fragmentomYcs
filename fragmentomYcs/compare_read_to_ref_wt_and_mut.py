from __future__ import annotations


def compare_read_to_ref_wt_and_mut(
    read_seq: str,
    ref_seq_wt: str,
    ref_seq_mut: str,
    compare_len_wt: int,
    compare_len_mut: int,
) -> str:
    ref_seq_wt_sub = ref_seq_wt[:compare_len_wt]
    ref_seq_mut_sub = ref_seq_mut[:compare_len_mut]
    read_seq_wt_sub = read_seq[:compare_len_wt]
    read_seq_mut_sub = read_seq[:compare_len_mut]

    match_ref_wt = ref_seq_wt_sub == read_seq_wt_sub
    match_ref_mut = ref_seq_mut_sub == read_seq_mut_sub

    if match_ref_wt and match_ref_mut:
        return "AMB"
    if match_ref_wt:
        return "WT"
    if match_ref_mut:
        return "MUT"
    return "OTH"
