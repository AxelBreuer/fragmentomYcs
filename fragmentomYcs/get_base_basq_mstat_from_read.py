from __future__ import annotations

from typing import Any, Dict, Mapping, Optional

from fragmentomYcs.get_index_aligning_with_pos import get_index_aligning_with_pos
from fragmentomYcs.get_mutation_status_of_read import get_mutation_status_of_read


def get_base_basq_mstat_from_read(
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    read_stats: Mapping[str, Any],
    fasta_seq: Mapping[str, Any],
) -> Dict[str, Optional[str]]:
    # get index in the read sequence aligning with pos
    read_index_at_pos = get_index_aligning_with_pos(pos, read_stats)

    if read_index_at_pos in (-1.0, -0.5):
        # in case the read does not cover the position of interest
        return {"base": None, "basq": None, "mstat": None}

    if read_index_at_pos == -2.0:
        # in case the read contains a deletion at the position of interest
        # then the read contains another event than the one we are looking
        # for
        mstat = "OTH"
    else:
        # get mutation status of the read
        if len(ref) == len(alt):
            # SNV or MNV
            mstat_small = get_mutation_status_of_read(
                chr_,
                pos,
                ref,
                alt,
                read_stats,
                int(read_index_at_pos),
                fasta_seq,
                n_match_base_before=0,
                n_match_base_after=0,
            )
            mstat_large = get_mutation_status_of_read(
                chr_,
                pos,
                ref,
                alt,
                read_stats,
                int(read_index_at_pos),
                fasta_seq,
                n_match_base_before=1,
                n_match_base_after=1,
            )

            if mstat_small == "MUT":
                if mstat_large == "OTH":
                    mstat = "MUT but potentially OTH"
                elif mstat_large == "MUT":
                    mstat = "MUT"
                else:
                    raise ValueError(
                        "If mutation status on the alt sequence is MUT, then extended "
                        f"status cannot be '{mstat_large}' for {chr_}:{pos}:{ref}>{alt}"
                    )
            else:
                mstat = mstat_small
        else:
            # By VCF convention, the ref and alt sequences include one base
            # before the actual mutated sequence. Additionally, we want to
            # systematically include one base after the actual mutated
            # sequences.
            mstat = get_mutation_status_of_read(
                chr_,
                pos,
                ref,
                alt,
                read_stats,
                int(read_index_at_pos),
                fasta_seq,
                n_match_base_before=1,
                n_match_base_after=1,
            )

    # get bases to report. The idea is to report the bases aligning with
    # pos, [inserted bases], pos+1, ..., pos+len(alt)
    # and stop reporting bases once we reach len(alt)

    # initialize
    if read_index_at_pos == -2.0:
        base = "-"
        basq = ""
    else:
        idx0 = int(read_index_at_pos) - 1
        seq = str(read_stats["SEQ"])
        qual = str(read_stats["QUAL"])
        base = seq[idx0 : idx0 + 1]
        basq = qual[idx0 : idx0 + 1]

    # iterate
    shift = 0
    while len(base) < len(alt):
        shift += 1
        info_at_pos = get_base_basq_from_read_at_pos(pos + shift, pos + shift - 1, read_stats)

        max_new = min(len(alt) - len(base), len(info_at_pos["base"]))
        base += info_at_pos["base"][:max_new]
        basq += info_at_pos["basq"][:max_new]

        if info_at_pos["base"] == "*":
            break

    return {"base": base, "basq": basq, "mstat": mstat}


def get_base_basq_from_read_at_pos(
    pos_cur: int,
    pos_pre: int,
    read_stats: Mapping[str, Any],
) -> Dict[str, str]:
    read_index_at_pos_cur = get_index_aligning_with_pos(pos_cur, read_stats)
    read_index_at_pos_pre = get_index_aligning_with_pos(pos_pre, read_stats)

    seq = str(read_stats["SEQ"])
    qual = str(read_stats["QUAL"])

    if read_index_at_pos_cur == -1.0:
        return {"base": "*", "basq": "*"}

    if read_index_at_pos_cur == -2.0:
        return {"base": "-", "basq": ""}

    if read_index_at_pos_cur == -0.5:
        start_idx = int(read_index_at_pos_pre)
        return {"base": seq[start_idx:], "basq": qual[start_idx:]}

    start_idx = int(read_index_at_pos_pre)
    end_idx = int(read_index_at_pos_cur)
    return {"base": seq[start_idx:end_idx], "basq": qual[start_idx:end_idx]}
