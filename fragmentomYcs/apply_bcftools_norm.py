"""Dispatcher for VCF variant left-alignment.

Reads ``[bcftools] backend`` from ``fragmentomYcs.cfg`` to select:

- ``auto`` *(default)* — Windows → pybcftools; Linux/macOS → bcftools binary
  if on PATH, otherwise pybcftools.
- ``pybcftools`` — pure-Python, no external tool, cross-platform.
- ``bcftools`` — system binary (Linux / WSL).

Edit ``fragmentomYcs.cfg`` to switch::

    [bcftools]
    backend = pybcftools   # or: bcftools, auto
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from typing import Any, Dict, Optional

import pyfaidx

from fragmentomYcs.config import get_bcftools_backend
from fragmentomYcs.pybcftools import left_align_and_trim


# ---------------------------------------------------------------------------
# Backend detection
# ---------------------------------------------------------------------------

def get_backend() -> str:
    """Return ``"pybcftools"`` or ``"bcftools"`` based on ``fragmentomYcs.cfg``."""
    cfg_val = get_bcftools_backend()  # "auto", "pybcftools", or "bcftools"

    if cfg_val == "pybcftools":
        return "pybcftools"
    if cfg_val == "bcftools":
        return "bcftools"

    # auto: Windows always uses pybcftools; elsewhere prefer binary if available
    if sys.platform == "win32":
        return "pybcftools"
    return "bcftools" if shutil.which("bcftools") is not None else "pybcftools"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def apply_bcftools_norm(
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    fasta: pyfaidx.Fasta,
    fasta_path: Optional[str] = None,
    verbose: bool = False,
) -> Optional[Dict[str, Any]]:
    """Normalise a biallelic VCF variant (left-align + minimal representation).

    Dispatches to the backend selected by ``FRAGMENTOMYCS_BCFTOOLS_BACKEND``
    (see module docstring).

    Parameters
    ----------
    chr_ : str
    pos  : int  1-based genomic position
    ref  : str  VCF ref allele (anchor base included for indels)
    alt  : str  VCF alt allele
    fasta : pyfaidx.Fasta
        Open FASTA handle — used by the ``pybcftools`` backend.
    fasta_path : str, optional
        Filesystem path to the FASTA file — required by the ``bcftools``
        binary backend.  Ignored when ``pybcftools`` is used.
    verbose : bool

    Returns
    -------
    dict with keys ``chr``, ``pos``, ``ref``, ``alt``, or ``None`` on error.
    """
    if ref == alt:
        return {"chr": chr_, "pos": pos, "ref": ref, "alt": alt}

    backend = get_backend()

    # Fall back to pybcftools if the binary backend cannot be satisfied
    if backend == "bcftools":
        reason = (
            "fasta_path not provided" if not fasta_path
            else "bcftools not on PATH" if shutil.which("bcftools") is None
            else ""
        )
        if reason:
            warnings.warn(f"{reason}; falling back to pybcftools.", stacklevel=2)
            backend = "pybcftools"

    if verbose:
        print(f"  [{backend}] Left-aligning {chr_}:{pos}:{ref}>{alt}")

    if backend == "bcftools":
        return _apply_bcftools_binary(chr_, pos, ref, alt, fasta_path, verbose)  # type: ignore[arg-type]

    # pybcftools backend
    try:
        new_pos, new_ref, new_alt = left_align_and_trim(chr_, pos, ref, alt, fasta)
        return {"chr": chr_, "pos": new_pos, "ref": new_ref, "alt": new_alt}
    except Exception as exc:
        warnings.warn(f"pybcftools failed for {chr_}:{pos}:{ref}>{alt}: {exc}", stacklevel=2)
        return None


# ---------------------------------------------------------------------------
# bcftools binary backend
# ---------------------------------------------------------------------------

def _apply_bcftools_binary(
    chr_: str,
    pos: int,
    ref: str,
    alt: str,
    fasta_path: str,
    verbose: bool,
) -> Optional[Dict[str, Any]]:
    """Call the system ``bcftools norm`` binary and parse the normalised result."""
    tmp_in = _write_temp_vcf(chr_, pos, ref, alt)
    tmp_out_fd, tmp_out = tempfile.mkstemp(suffix=".vcf")
    os.close(tmp_out_fd)

    try:
        cmd = [
            "bcftools", "norm",
            "-m", "+both",
            "-d", "exact",
            "--check", "REF,ALT",
            "-f", fasta_path,
            "-o", tmp_out,
            tmp_in,
        ]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode != 0:
            warnings.warn(
                f"bcftools norm failed for {chr_}:{pos}:{ref}>{alt}: "
                f"{result.stderr.decode().strip()}",
                stacklevel=3,
            )
            return None

        if not os.path.isfile(tmp_out) or os.path.getsize(tmp_out) == 0:
            warnings.warn(
                f"bcftools norm produced an empty VCF for {chr_}:{pos}:{ref}>{alt}.",
                stacklevel=3,
            )
            return None

        rows = _parse_vcf(tmp_out)
        if not rows:
            warnings.warn(
                f"bcftools norm returned no variants for {chr_}:{pos}:{ref}>{alt}.",
                stacklevel=3,
            )
            return None

        r = rows[0]
        return {
            "chr": str(r["CHROM"]),
            "pos": int(r["POS"]),
            "ref": str(r["REF"]),
            "alt": str(r["ALT"]),
        }
    finally:
        for path in (tmp_in, tmp_out):
            try:
                os.remove(path)
            except OSError:
                pass


def _write_temp_vcf(chr_: str, pos: int, ref: str, alt: str) -> str:
    """Write a minimal single-variant VCF to a temp file; return its path."""
    header = (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={chr_}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    record = f"{chr_}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n"
    fd, path = tempfile.mkstemp(suffix=".vcf")
    with os.fdopen(fd, "w") as fh:
        fh.write(header + record)
    return path


def _parse_vcf(vcf_path: str) -> list:
    """Parse CHROM/POS/REF/ALT from a VCF file (skips ``#`` header lines)."""
    rows = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            rows.append({
                "CHROM": parts[0],
                "POS": int(parts[1]),
                "REF": parts[3],
                "ALT": parts[4],
            })
    return rows
