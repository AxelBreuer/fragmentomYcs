from __future__ import annotations

import re
from typing import Dict


def get_fragment_bases_5p_3p_softclip(cigar_5p: str, cigar_3p: str) -> Dict[str, int]:
    m_5p = re.match(r"^(\d+)S", cigar_5p)
    m_3p = re.search(r"(\d+)S$", cigar_3p)

    return {
        "nb_softclip_5p": int(m_5p.group(1)) if m_5p else 0,
        "nb_softclip_3p": int(m_3p.group(1)) if m_3p else 0,
    }
