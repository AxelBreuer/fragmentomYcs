from __future__ import annotations

from typing import Any, Mapping, Sequence

from fragmentomYcs.utils import is_na, sam_flag_bits


def process_fragment_reads_qc(df_fragment_reads: Sequence[Mapping[str, Any]], chr_: str) -> str:
    qc_messages = []
    num_reads = len(df_fragment_reads)

    # --- CASE A: Abnormal Read Count (Not 2) ---
    if num_reads != 2:
        # Detailed diagnosis if we have exactly 1 read
        if num_reads == 1:
            read = df_fragment_reads[0]
            msg = "Fragment has 1 read"

            # Normalize mate chromosome info
            mate_chr = read.get("RNEXT")
            if not is_na(mate_chr) and mate_chr == "=":
                mate_chr = read.get("RNAME")

            # Diagnosis 1: Mate is unmapped (RNEXT is *)
            if is_na(mate_chr) or mate_chr == "*":
                msg = msg + " - Mate is unmapped"
            # Diagnosis 2: Translocation (Mate on diff chr, filtered out by read_bam)
            elif mate_chr != read.get("RNAME"):
                msg = msg + f" - Mate maps to a different chromosome: {mate_chr}"
            # Diagnosis 3: Far Away (Mate on same chr but outside window, filtered out)
            else:
                tlen = read.get("TLEN")
                dist = "NA" if is_na(tlen) else str(abs(int(tlen)))
                msg = msg + f" - Mate maps outside loaded region. TLEN: {dist}"
        else:
            msg = f"Fragment has {num_reads} read(s)"

        qc_messages.append(msg)

    # --- CASE B: Normal Read Count (2 Reads) - Check Consistency ---
    if num_reads == 2:
        # Use first read to check pairing info consistency
        read = df_fragment_reads[0]

        # Normalize mate chromosome
        mate_chr = read.get("RNEXT")
        if not is_na(mate_chr) and mate_chr == "=":
            mate_chr = read.get("RNAME")

        # Test 1: Wrong Chromosome (Are the loaded reads actually on the target chr?)
        if not all(x.get("RNAME") == chr_ for x in df_fragment_reads):
            qc_messages.append(f"Read(s) found on a chromosome other than {chr_}")

        # Test 2: Missing Mate Info (RNEXT is *)
        if is_na(read.get("RNEXT")) or read.get("RNEXT") == "*":
            qc_messages.append("Mate chromosome info is not available")
        # Test 3: Internal Consistency (Do reads claim to be on different chromosomes?)
        elif mate_chr != read.get("RNAME"):
            qc_messages.append(f"- Mate maps to a different chromosome: {mate_chr}")

        # Test 4: Mapping Status
        if any(is_na(x.get("POS")) or int(x.get("POS", 0)) == 0 for x in df_fragment_reads):
            qc_messages.append("One or both reads are unmapped (POS=0 or NA)")

        # Test 5: Valid Pair (Must contain exactly one R1 and one R2)
        # Test 6: Strand Orientation (Must have one forward and one reverse read)
        flag_matrix = [sam_flag_bits(int(x["FLAG"])) for x in df_fragment_reads]

        if sum(1 for bits in flag_matrix if bits["isFirstMateRead"]) != 1:
            qc_messages.append("Fragment is not a valid R1/R2 pair")

        if sum(1 for bits in flag_matrix if bits["isMinusStrand"]) != 1:
            qc_messages.append("Badly oriented reads")

    # --- CASE C: Strand Orientation for single-read fragments ---
    # When the mate was not loaded, we can still infer orientation from the
    # FLAG of the available read: it encodes both its own strand (isMinusStrand)
    # and its mate's strand (isMateMinusStrand).
    if num_reads == 1:
        bits = sam_flag_bits(int(df_fragment_reads[0]["FLAG"]))
        if bits["isMinusStrand"] == bits["isMateMinusStrand"]:
            qc_messages.append("Badly oriented reads")

    return " & ".join(qc_messages)
