from __future__ import annotations

import re
from typing import Any, Mapping

from fragmentomYcs.compare_read_to_ref_wt_and_mut import compare_read_to_ref_wt_and_mut
from fragmentomYcs.search_for_indel_in_cigar import search_for_indel_in_cigar
from fragmentomYcs.utils import get_number_of_common_first_char
from fragmentomYcs.utils_fasta import get_seq_from_fasta


def get_mutation_status_of_read(
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    read_stats: Mapping[str, Any],
    read_index_at_pos: int,
    fasta_seq: Mapping[str, Any],
    n_match_base_before: int = 1,
    n_match_base_after: int = 1,
) -> str:
    ref_len = len(ref)
    alt_len = len(alt)
    read_seq = str(read_stats["SEQ"])
    read_seq_len = len(read_seq)

    # first identify the minimum number of bases that should be considered in
    # the comparison
    if ref_len == alt_len:
        # SNV or MNV

        # If the read sequence does not allow to compare with base(s) before
        # the alteration, or with base(s) after the alteration, then we won't
        # add these bases into the comparison
        if read_index_at_pos == 1:
            n_match_base_before = 0
        if (read_index_at_pos + alt_len - 1) >= read_seq_len:
            n_match_base_after = 0

        # we fetch one base before, the size of the allele, and one base after
        fetch_len_ref = n_match_base_before + alt_len + n_match_base_after
        fetch_start_ref = pos - n_match_base_before
        fetch_end_ref = fetch_start_ref + fetch_len_ref - 1

        # get sequences of wild-type and mutated ref
        ref_seq_wt = get_seq_from_fasta(chr_, fetch_start_ref, fetch_end_ref, fasta_seq=fasta_seq)
        ref_seq_mut = (
            ref_seq_wt[:n_match_base_before]
            + alt
            + ref_seq_wt[len(ref_seq_wt) - n_match_base_after :]
        )

        # if we need more bases than available, cut the sequences
        fetch_start_read = read_index_at_pos - n_match_base_before
        start_idx = max(fetch_start_read - 1, 0)
        fetch_len_read = min(fetch_len_ref, read_seq_len - start_idx)
        fetch_end_read = start_idx + fetch_len_read
        read_subseq = read_seq[start_idx:fetch_end_read]

        # run comparison on maximum size possible
        compare_len_wt = n_match_base_before + alt_len + n_match_base_after
        compare_len_mut = n_match_base_before + alt_len + n_match_base_after

        # if we cannot cover completely the allele with the read, we may have
        # an ambiguity
        incomplete_comparison_mut = fetch_len_read < compare_len_mut
        compare_len_wt = min(compare_len_wt, fetch_len_read)
        compare_len_mut = min(compare_len_mut, fetch_len_read)

        # determine the read mutation status
        status = compare_read_to_ref_wt_and_mut(
            read_subseq,
            ref_seq_wt,
            ref_seq_mut,
            compare_len_wt,
            compare_len_mut,
        )

        if status == "MUT" and incomplete_comparison_mut:
            # An incomplete comparison compatible with the mutated allele
            # will be labelled 'AMB'
            return "AMB"
        return status

    # INS or DEL
    # identify mutated sequence
    inserted_seq = alt[1:] if alt_len > ref_len else ""
    deleted_seq = ref[1:] if ref_len > alt_len else ""
    motif = inserted_seq if inserted_seq else deleted_seq
    motif_len = len(motif)

    # Fetch the maximum number of bases that may be included in the
    # comparison
    fetch_len_ref = n_match_base_before - 1 + read_seq_len - (read_index_at_pos - 1) + motif_len + n_match_base_after
    fetch_start_ref = pos - (n_match_base_before - 1)
    fetch_end_ref = fetch_start_ref + fetch_len_ref - 1

    ref_seq_wt = get_seq_from_fasta(chr_, fetch_start_ref, fetch_end_ref, fasta_seq=fasta_seq)

    # if the sequence contains the searched motif at least once
    regex_motif = "^" + (alt if alt_len > ref_len else ref)
    if re.search(regex_motif, ref_seq_wt[n_match_base_before - 1 :]):
        repeat_count = 1
        ref_seq_after_one_motif = ref_seq_wt[n_match_base_before + motif_len :]

        # identify the number of repeats within the region of the maximum
        # comparison size
        while ref_seq_after_one_motif[(repeat_count - 1) * motif_len : repeat_count * motif_len] == motif:
            repeat_count += 1

        start_after_last = (repeat_count - 1) * motif_len + (n_match_base_after - 1)
        end_after_last = repeat_count * motif_len + (n_match_base_after - 1)
        ref_seq_after_last_motif = ref_seq_after_one_motif[start_after_last:end_after_last]
        n_bases_shared_with_motif = get_number_of_common_first_char(ref_seq_after_last_motif, motif)
    else:
        repeat_count = 0
        n_bases_shared_with_motif = 0

    # Compute the size of the comparisons to wild-type ref and mutated ref
    if alt_len > ref_len:
        indel_type = "I"
        compare_len_wt = n_match_base_before + repeat_count * motif_len + n_bases_shared_with_motif + n_match_base_after
        compare_len_mut = n_match_base_before + (repeat_count + 1) * motif_len + n_bases_shared_with_motif + n_match_base_after
    else:
        indel_type = "D"
        compare_len_wt = n_match_base_before + (repeat_count - 1) * motif_len + n_bases_shared_with_motif + n_match_base_after
        compare_len_mut = n_match_base_before + (repeat_count - 1) * motif_len + n_bases_shared_with_motif + n_match_base_after

    # Build mutated sequence by inserting or deleting the correct bases
    ref_seq_mut = ref_seq_wt[:n_match_base_before] + alt[n_match_base_before:] + ref_seq_wt[ref_len:]

    # Fetch maximum read sequence available starting from the position
    # needed to cover n_match_base_before
    fetch_start_read = read_index_at_pos - (n_match_base_before - 1)
    start_idx = max(fetch_start_read - 1, 0)
    fetch_len_read = read_seq_len - start_idx
    read_subseq = read_seq[start_idx : start_idx + fetch_len_read]

    # if we cannot cover completely the allele with the read, we may have
    # an ambiguity
    incomplete_comparison_mut = fetch_len_read < compare_len_mut
    compare_len_wt = min(compare_len_wt, fetch_len_read)
    compare_len_mut = min(compare_len_mut, fetch_len_read)

    # determine the read mutation status
    status = compare_read_to_ref_wt_and_mut(
        read_subseq,
        ref_seq_wt,
        ref_seq_mut,
        compare_len_wt,
        compare_len_mut,
    )

    indel_found_in_cigar, other_found_in_cigar = search_for_indel_in_cigar(
        pos,
        ref,
        alt,
        read_stats,
        indel_type,
    )

    if not incomplete_comparison_mut:
        # -------------------- complete_comparison --------------------
        if indel_found_in_cigar:
            return {
                "MUT": "MUT",
                "WT": "MUT by CIGAR but potentially WT",
                "AMB": "IMPOSSIBLE",
                "OTH": "MUT by CIGAR but potentially OTH",
            }[status]

        if other_found_in_cigar:
            return {
                "MUT": "OTH by CIGAR but potentially MUT",
                "WT": "OTH by CIGAR but potentially WT",
                "AMB": "IMPOSSIBLE",
                "OTH": "OTH",
            }[status]

        return {
            "MUT": "WT by CIGAR but potentially MUT",
            "WT": "WT",
            "AMB": "IMPOSSIBLE",
            "OTH": "OTH",
        }[status]

    # -------------------- incomplete_comparison --------------------
    if indel_found_in_cigar:
        return {
            "MUT": "MUT by CIGAR but AMB",
            "WT": "MUT by CIGAR but potentially WT",
            "AMB": "MUT by CIGAR but AMB",
            "OTH": "MUT by CIGAR but potentially OTH",
        }[status]

    if other_found_in_cigar:
        return {
            "MUT": "OTH by CIGAR but potentially MUT",
            "WT": "OTH by CIGAR but potentially WT",
            "AMB": "OTH by CIGAR but AMB",
            "OTH": "OTH",
        }[status]

    return {
        "MUT": "WT by CIGAR but potentially MUT",
        "WT": "WT",
        "AMB": "AMB",
        "OTH": "OTH",
    }[status]
