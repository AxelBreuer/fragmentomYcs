from __future__ import annotations

import re
import warnings
from typing import List, Dict, Any

import pandas as pd


def read_mut(mut: str) -> List[Dict[str, Any]]:
    """Read variant information from multiple sources.

    Acts as a dispatcher to read variant information from a file (VCF or TSV)
    or a "chr:pos:ref:alt" string.  All multi-allelic sites are split into
    separate bi-allelic rows.

    Mirrors read_mut() / its sub-routines in R/read_mut.R.

    Parameters
    ----------
    mut:
        One of:
        * Path to a TSV file with columns CHROM, POS, REF, ALT
          (plain or ``.gz`` compressed).
        * Path to a VCF file (plain or ``.gz`` / bgzf compressed).
        * A ``"chr:pos:ref:alt"`` formatted string.

    Returns
    -------
    list[dict]
        Each dict has keys ``CHROM``, ``POS`` (int), ``REF``, ``ALT``.
    """
    if re.search(r"\.tsv(\.gz)?$", mut):
        df = _read_tsv_input(mut)
    elif re.search(r"\.vcf(\.gz)?$", mut):
        df = _read_vcf_input(mut)
    elif re.search(r"^[^:]+:[^:]+:[^:]+:[^:]+$", mut):
        df = _read_str_input(mut)
    else:
        raise ValueError(
            f"The parameter 'mut' ('{mut}') is not in the expected format "
            "(.tsv, .vcf, chr:pos:ref:alt)."
        )

    return df.to_dict(orient="records")


# ---------------------------------------------------------------------------
# Top-level dispatched parsers
# ---------------------------------------------------------------------------

def _read_vcf_input(vcf_file: str) -> pd.DataFrame:
    """Read VCF and expand multi-allelic records."""
    df = _parser_vcf(vcf_file)
    return _expand_multiallelics(df)


def _read_tsv_input(tsv_file: str) -> pd.DataFrame:
    """Read TSV and expand multi-allelic records."""
    df = _parser_tsv(tsv_file)
    return _expand_multiallelics(df)


def _read_str_input(mut: str) -> pd.DataFrame:
    """Parse a chr:pos:ref:alt string and expand multi-allelic records."""
    df = _parser_str(mut)
    return _expand_multiallelics(df)


# ---------------------------------------------------------------------------
# Low-level parsers
# ---------------------------------------------------------------------------

def _parser_vcf(vcf_file: str) -> pd.DataFrame:
    """Parse CHROM/POS/REF/ALT from a VCF file using cyvcf2.

    Falls back to pysam.VariantFile if cyvcf2 is not installed.
    """
    try:
        import cyvcf2  # preferred: fast, handles bgzf and plain gzip
        return _parser_vcf_cyvcf2(vcf_file, cyvcf2)
    except ImportError:
        pass

    try:
        import pysam  # fallback
        return _parser_vcf_pysam(vcf_file, pysam)
    except ImportError:
        raise ImportError(
            "A VCF parser is required. Install cyvcf2 (recommended) or pysam:\n"
            "  pip install cyvcf2\n"
            "  conda install -c bioconda cyvcf2"
        )


def _parser_vcf_cyvcf2(vcf_file: str, cyvcf2) -> pd.DataFrame:
    rows = []
    vcf = cyvcf2.VCF(vcf_file)
    for variant in vcf:
        chrom = variant.CHROM
        pos = variant.POS  # 1-based
        ref = variant.REF
        alt = ",".join(variant.ALT) if variant.ALT else "-"
        rows.append({"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt})
    vcf.close()

    if not rows:
        raise ValueError("The VCF file does not contain any mutation data.")

    return pd.DataFrame(rows)


def _parser_vcf_pysam(vcf_file: str, pysam) -> pd.DataFrame:
    rows = []
    vcf = pysam.VariantFile(vcf_file)
    for rec in vcf:
        chrom = rec.chrom
        pos = rec.pos + 1  # pysam is 0-based
        ref = rec.ref
        alt = ",".join(rec.alts) if rec.alts else "-"
        rows.append({"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt})
    vcf.close()

    if not rows:
        raise ValueError("The VCF file does not contain any mutation data.")

    return pd.DataFrame(rows)


def _parser_tsv(tsv_file: str) -> pd.DataFrame:
    """Read a TSV file with columns CHROM, POS, REF, ALT."""
    try:
        df = pd.read_csv(
            tsv_file,
            sep="\t",
            usecols=["CHROM", "POS", "REF", "ALT"],
            dtype={"CHROM": str, "POS": str, "REF": str, "ALT": str},
            keep_default_na=False,
        )
    except Exception as e:
        raise ValueError("Failed to read the TSV file. Ensure it's properly formatted.") from e

    if df.shape[1] != 4:
        raise ValueError(
            f"The number of columns in the TSV file does not match the expected "
            f"columns CHROM POS REF ALT. Found columns: {', '.join(df.columns)}"
        )

    # Identify lines with non-integer values in POS
    invalid_pos_mask = ~df["POS"].str.match(r"^[0-9]+$") & df["POS"].notna()
    if invalid_pos_mask.any():
        invalid_rows = list(df.index[invalid_pos_mask] + 1)  # 1-based row numbers
        warnings.warn(
            f"The following rows have non-integer values in the POS column "
            f"and will be converted to NA: {', '.join(map(str, invalid_rows))}"
        )

    # Convert POS to integer (invalid become NaN then pd.NA)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce").astype("Int64")

    # Remove empty rows (all columns NA or empty string)
    is_empty = (df.isna() | (df == "")).all(axis=1)
    df = df[~is_empty].reset_index(drop=True)

    if df.empty:
        raise ValueError("The TSV file does not contain any mutation data.")

    return df


def _parser_str(mut: str) -> pd.DataFrame:
    """Parse a 'chr:pos:ref:alt' formatted string."""
    # If 'mut' ends with ":", append "-" to allow proper splitting
    if mut.endswith(":"):
        mut = mut + "-"

    parts = mut.split(":")

    if len(parts) != 4:
        raise ValueError(
            f"The parameter 'mut' ('{mut}') is not in the expected format "
            "(.tsv, .vcf, chr:pos:ref:alt)."
        )

    chr_ = parts[0]

    if not re.match(r"^[0-9]+$", parts[1]):
        warnings.warn(
            f"Position value '{parts[1]}' is not a valid integer and will be treated as NA."
        )
        pos = pd.NA
    else:
        pos = int(parts[1])

    ref = parts[2]
    alt = parts[3]

    df = pd.DataFrame([{"CHROM": chr_, "POS": pos, "REF": ref, "ALT": alt}])

    if df.empty:
        raise ValueError("The mutation string does not contain any mutation data.")

    return df


# ---------------------------------------------------------------------------
# Multi-allelic expansion
# ---------------------------------------------------------------------------

def _expand_multiallelics(df: pd.DataFrame) -> pd.DataFrame:
    """Expand rows where ALT contains comma-separated alleles.

    Mirrors expand_multiallelics() in R/read_mut.R.
    """
    # Identify and remove rows where REF itself is multi-allelic
    is_multi_ref = df["REF"].str.contains(",", na=False)
    if is_multi_ref.any():
        warnings.warn("REF can not be multiallelic. These variants have been removed.")
        df = df[~is_multi_ref].reset_index(drop=True)

    if df.empty:
        raise ValueError(
            "After reading and filtering, the mutation information is empty. "
            "No valid mutation data found."
        )

    # Replace NA or empty string with "-"
    df["ALT"] = df["ALT"].fillna("-").replace("", "-")

    # Append a "-" if ALT ends with a comma (e.g. "A," becomes "A,-")
    df["ALT"] = df["ALT"].str.replace(r",$", ",-", regex=True)

    # Split on comma and explode (equivalent to tidyr::separate_rows)
    df = df.assign(ALT=df["ALT"].str.split(",")).explode("ALT").reset_index(drop=True)

    return df
