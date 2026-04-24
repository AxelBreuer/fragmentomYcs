from __future__ import annotations

from typing import Any, Dict, Mapping, Optional

from fragmentomYcs.get_base_basq_mstat_from_read import get_base_basq_mstat_from_read
from fragmentomYcs.get_fragment_bases_5p_3p import get_fragment_bases_5p_3p
from fragmentomYcs.get_fragment_bases_5p_3p_softclip import get_fragment_bases_5p_3p_softclip
from fragmentomYcs.get_fragment_size import get_fragment_size
from fragmentomYcs.get_mutation_status_of_fragment import get_mutation_status_of_fragment
from fragmentomYcs.process_fragment_reads_QC import process_fragment_reads_qc
from fragmentomYcs.remove_softclip import remove_softclip as remove_softclip_fn
from fragmentomYcs.utils import as_records, is_na, parse_cigar, return_fail_qc_fragment, sam_flag_bits


def extract_fragment_features(
    df_sam: Any,
    fragment_name: str,
    sample_id: Optional[str],
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    report_bam_info: bool,
    report_softclip: bool,
    report_5p_3p_bases_fragment: int,
    remove_softclip: bool,
    fasta_seq: Mapping[str, Any],
    input_mutation_info: str,
) -> Dict[str, Any]:
    # ----- Read QCs -----
    # Select reads originating from the fragment of interest
    mask_frag = [row for row in as_records(df_sam) if row.get("QNAME") == fragment_name]
    df_fragment_reads = mask_frag

    # Sanity check fragments
    fragment_qc = process_fragment_reads_qc(df_fragment_reads, chr_)

    # If the fragment fails QC, return a dict with Nones
    if fragment_qc != "":
        return return_fail_qc_fragment(
            fragment_qc,
            sample_id,
            chr_,
            pos,
            ref,
            alt,
            input_mutation_info,
            fragment_name,
        )

    # ----- Read preprocessing -----
    # Define 3' and 5' reads
    # Get FLAG attributes for both reads in the fragment.
    flag_matrix = [sam_flag_bits(int(r["FLAG"])) for r in df_fragment_reads]

    # Identify the row index of the 5p read (forward strand, where
    # isMinusStrand is 0).
    idx_5p = next(i for i, bits in enumerate(flag_matrix) if not bits["isMinusStrand"])

    # The 3p read is the other one (reverse strand, where isMinusStrand is 1).
    idx_3p = next(i for i, bits in enumerate(flag_matrix) if bits["isMinusStrand"])

    # Get read bam info for read 5' and read 3'
    read_stats_5p = get_read_stats(df_fragment_reads[idx_5p])
    read_stats_3p = get_read_stats(df_fragment_reads[idx_3p])

    # If remove_softclip = True, remove
    # softclip for analysis
    if remove_softclip:
        # Keep in memory the input information
        input_read_stats_5p = dict(read_stats_5p)
        input_read_stats_3p = dict(read_stats_3p)

        # Remove softclip
        read_5p_info_without_softclip = remove_softclip_fn(read_stats_5p)
        read_3p_info_without_softclip = remove_softclip_fn(read_stats_3p)

        if read_5p_info_without_softclip["SEQ"] == "" or read_5p_info_without_softclip["CIGAR"] == "":
            return return_fail_qc_fragment(
                "Invalid read after softclip trimming",
                sample_id,
                chr_,
                pos,
                ref,
                alt,
                input_mutation_info,
                fragment_name,
            )

        if read_3p_info_without_softclip["SEQ"] == "" or read_3p_info_without_softclip["CIGAR"] == "":
            return return_fail_qc_fragment(
                "Invalid read after softclip trimming",
                sample_id,
                chr_,
                pos,
                ref,
                alt,
                input_mutation_info,
                fragment_name,
            )

        read_stats_5p.update(read_5p_info_without_softclip)
        read_stats_3p.update(read_3p_info_without_softclip)
    else:
        input_read_stats_5p = read_stats_5p
        input_read_stats_3p = read_stats_3p

    # ----- Fragmentomic features extraction -----
    # Get read sequence, read base qualities,
    # and read mutation status
    read_info_5p = get_base_basq_mstat_from_read(chr_, pos, ref, alt, read_stats_5p, fasta_seq)
    read_info_3p = get_base_basq_mstat_from_read(chr_, pos, ref, alt, read_stats_3p, fasta_seq)

    # Compute fragment size
    fragment_size = get_fragment_size(read_stats_5p, read_stats_3p)

    # Define fragment status
    fstats = get_mutation_status_of_fragment(
        mstat_5p=read_info_5p["mstat"],
        mstat_3p=read_info_3p["mstat"],
    )

    # Compute Position_3p -> last aligned position of the fragment
    Position_3p = end_on_reference(read_stats_3p["POS"], read_stats_3p["CIGAR"])

    # Build an adaptative output
    output_read_stats_5p = input_read_stats_5p if remove_softclip else read_stats_5p
    output_read_stats_3p = input_read_stats_3p if remove_softclip else read_stats_3p

    final_row_fragment: Dict[str, Any] = {
        "Sample_Id": sample_id,
        "Chromosome": chr_,
        "Position": int(pos),
        "Ref": ref,
        "Alt": alt,
        "Input_Mutation": input_mutation_info,
        "Fragment_Id": fragment_name,
        "Fragment_QC": "OK",
        "Fragment_Status_Simple": fstats["Simple"],
        "Fragment_Status_Detail": fstats["Detail"],
        "Fragment_Size": int(fragment_size),
        "Read_5p_Status": read_info_5p["mstat"],
        "Read_3p_Status": read_info_3p["mstat"],
        "BASE_5p": read_info_5p["base"],
        "BASE_3p": read_info_3p["base"],
        "BASQ_5p": read_info_5p["basq"],
        "BASQ_3p": read_info_3p["basq"],
        "Position_5p": int(output_read_stats_5p["POS"]),
        "Position_3p": int(Position_3p) if Position_3p is not None else None,
    }

    if report_bam_info:
        final_row_fragment["POS_5p"] = int(output_read_stats_5p["POS"])
        final_row_fragment["POS_3p"] = int(output_read_stats_3p["POS"])
        final_row_fragment["FLAG_5p"] = int(output_read_stats_5p["FLAG"])
        final_row_fragment["FLAG_3p"] = int(output_read_stats_3p["FLAG"])
        final_row_fragment["MAPQ_5p"] = int(output_read_stats_5p["MAPQ"])
        final_row_fragment["MAPQ_3p"] = int(output_read_stats_3p["MAPQ"])
        final_row_fragment["CIGAR_5p"] = str(output_read_stats_5p["CIGAR"])
        final_row_fragment["CIGAR_3p"] = str(output_read_stats_3p["CIGAR"])

        tlen = read_stats_5p.get("TLEN")
        final_row_fragment["TLEN"] = None if is_na(tlen) else int(abs(int(tlen)))

    if int(report_5p_3p_bases_fragment) != 0:
        fragment_bases_5p_3p = get_fragment_bases_5p_3p(
            int(report_5p_3p_bases_fragment),
            read_stats_5p["SEQ"],
            read_stats_3p["SEQ"],
            read_stats_5p["QUAL"],
            read_stats_3p["QUAL"],
        )

        final_row_fragment["Fragment_Bases_5p"] = fragment_bases_5p_3p["fragment_bases_5p"]
        final_row_fragment["Fragment_Bases_3p"] = fragment_bases_5p_3p["fragment_bases_3p"]
        final_row_fragment["Fragment_Basqs_5p"] = fragment_bases_5p_3p["fragment_basqs_5p"]
        final_row_fragment["Fragment_Basqs_3p"] = fragment_bases_5p_3p["fragment_basqs_3p"]

    # Define number of soft clipped bases in 5'
    # If report_softclip is True, Call
    # the function get_fragment_bases_5p_3p_softclip
    if report_softclip:
        fragment_bases_softclip_5p_3p = get_fragment_bases_5p_3p_softclip(
            read_stats_5p["CIGAR"],
            read_stats_3p["CIGAR"],
        )

        final_row_fragment["Nb_Fragment_Bases_Softclip_5p"] = int(
            fragment_bases_softclip_5p_3p["nb_softclip_5p"]
        )
        final_row_fragment["Nb_Fragment_Bases_Softclip_3p"] = int(
            fragment_bases_softclip_5p_3p["nb_softclip_3p"]
        )

    # Add VAF column
    final_row_fragment["VAF"] = None

    return final_row_fragment


def get_read_stats(df_read: Mapping[str, Any]) -> Dict[str, Any]:
    seq = str(df_read["SEQ"])
    return {
        "FLAG": int(df_read["FLAG"]),
        "MAPQ": int(df_read["MAPQ"]),
        "TLEN": df_read.get("TLEN"),
        "CIGAR": str(df_read["CIGAR"]),
        "POS": int(df_read["POS"]),
        "SEQ": seq,
        "QUAL": str(df_read["QUAL"]),
        "read_length": len(seq),
    }


def end_on_reference(pos: int | None, cigar: str | None) -> int | None:
    if is_na(pos) or is_na(cigar) or cigar == "*":
        return None

    ops = parse_cigar(str(cigar))
    ref_len = sum(length for length, op in ops if op in ("M", "=", "X", "D", "N"))
    return int(pos) + int(ref_len) - 1
