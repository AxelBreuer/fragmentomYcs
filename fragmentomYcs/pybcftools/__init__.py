"""Pure-Python bcftools replacement.

Provides left-alignment of VCF variants without requiring the ``bcftools``
binary.  Suitable for Windows and any platform where ``bcftools`` is
unavailable.

Usage::

    from python.pybcftools import left_align_and_trim

    pos, ref, alt = left_align_and_trim("chr7", 55241864, "ATCG", "A", fasta)
"""

from fragmentomYcs.pybcftools.norm import left_align_and_trim

__all__ = ["left_align_and_trim"]
