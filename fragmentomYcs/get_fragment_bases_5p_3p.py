from __future__ import annotations

from typing import Dict


def get_fragment_bases_5p_3p(
    n_bases: int,
    seq_5p: str,
    seq_3p: str,
    qual_5p: str,
    qual_3p: str,
) -> Dict[str, str]:
    if n_bases > len(seq_5p):
        return {
            "fragment_bases_5p": seq_5p,
            "fragment_bases_3p": seq_3p,
            "fragment_basqs_5p": qual_5p,
            "fragment_basqs_3p": qual_3p,
        }

    return {
        "fragment_bases_5p": seq_5p[:n_bases],
        "fragment_bases_3p": seq_3p[-n_bases:],
        "fragment_basqs_5p": qual_5p[:n_bases],
        "fragment_basqs_3p": qual_3p[-n_bases:],
    }
