from __future__ import annotations

import re
from typing import Dict, List, Optional


def get_mutation_status_of_fragment(
    mstat_5p: Optional[str],
    mstat_3p: Optional[str],
) -> Dict[str, Optional[str]]:
    # Internal helper to clean status for logical comparisons (removes extra
    # detail)
    def clean_status(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        return re.sub(r"(:.*|\s.*)", "", s)

    # Function if a string contains "potentially X", extract X (WT/MUT/OTH/AMB);
    # else None
    def extract_potential_target(s: Optional[str]) -> Optional[str]:
        if s is None:
            return None
        m = re.search(r"potentially\s+(WT|MUT|OTH|AMB)", s, flags=re.IGNORECASE)
        return m.group(1).upper() if m else None

    # Function to combine original statuses for Fragment_Status_Detail
    # Handles sorting and unique values while preserving original string
    def combine_original_statuses(s1: Optional[str], s2: Optional[str]) -> Optional[str]:
        original_statuses = [s for s in (s1, s2) if s is not None]

        if len(original_statuses) == 0:
            return None
        if len(original_statuses) == 1:
            return original_statuses[0]
        if original_statuses[0] == original_statuses[1]:
            return original_statuses[0]

        pairs = [
            (clean_status(s1) or "", s1 or ""),
            (clean_status(s2) or "", s2 or ""),
        ]
        pairs.sort(key=lambda x: (x[0], x[1]))

        unique_sorted_originals: List[str] = []
        for _, original in pairs:
            if original and original not in unique_sorted_originals:
                unique_sorted_originals.append(original)

        return " & ".join(unique_sorted_originals)

    base_mstat_5p = clean_status(mstat_5p)
    base_mstat_3p = clean_status(mstat_3p)

    # Define flags for cleaner logic
    is_na1 = mstat_5p is None
    is_na2 = mstat_3p is None

    is_mut1 = base_mstat_5p == "MUT"
    is_mut2 = base_mstat_3p == "MUT"
    is_wt1 = base_mstat_5p == "WT"
    is_wt2 = base_mstat_3p == "WT"
    is_amb1 = base_mstat_5p == "AMB"
    is_amb2 = base_mstat_3p == "AMB"
    is_other_mut1 = base_mstat_5p == "OTH"
    is_other_mut2 = base_mstat_3p == "OTH"

    # --------------------------------------------------------------------------
    # Early resolution for discordant reads with "potentially X" on one side:
    # - If read A is "potentially X" and read B is X, trust B -> Simple = X
    # --------------------------------------------------------------------------
    if not is_na1 and not is_na2 and base_mstat_5p != base_mstat_3p:
        pot1 = extract_potential_target(mstat_5p)
        pot2 = extract_potential_target(mstat_3p)

        if pot1 is not None and base_mstat_3p == pot1:
            return {
                "Detail": combine_original_statuses(mstat_5p, mstat_3p),
                "Simple": base_mstat_3p,
            }

        if pot2 is not None and base_mstat_5p == pot2:
            return {
                "Detail": combine_original_statuses(mstat_5p, mstat_3p),
                "Simple": base_mstat_5p,
            }

    # --------------------------------------------------------------------------
    # Fragment status decision tree
    # --------------------------------------------------------------------------
    # 1. NA / NA
    if is_na1 and is_na2:
        return {"Detail": "ERR", "Simple": "ERR"}
    # 2. NA / MUT (and symmetrical)
    if (is_na1 and is_mut2) or (is_mut1 and is_na2):
        return {"Detail": mstat_3p if is_na1 else mstat_5p, "Simple": "MUT"}
    # 3. NA / WT (and symmetrical)
    if (is_na1 and is_wt2) or (is_wt1 and is_na2):
        return {"Detail": mstat_3p if is_na1 else mstat_5p, "Simple": "WT"}
    # 4. NA / AMB (and symmetrical)
    if (is_na1 and is_amb2) or (is_amb1 and is_na2):
        return {"Detail": mstat_3p if is_na1 else mstat_5p, "Simple": "N/I"}
    # 5. NA / OTH (and symmetrical)
    if (is_na1 and is_other_mut2) or (is_other_mut1 and is_na2):
        return {"Detail": mstat_3p if is_na1 else mstat_5p, "Simple": "OTH"}

    # 6. MUT / MUT
    if is_mut1 and is_mut2:
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "MUT"}
    # 7. MUT / WT (and symmetrical) - N/I
    if (is_mut1 and is_wt2) or (is_wt1 and is_mut2):
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "N/I"}
    # 8. MUT / AMB (and symmetrical) - Prioritize MUT
    if (is_mut1 and is_amb2) or (is_amb1 and is_mut2):
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "MUT"}
    if (is_mut1 and is_other_mut2) or (is_other_mut1 and is_mut2):
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "N/I"}

    if is_wt1 and is_wt2:
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "WT"}
    if (is_wt1 and is_amb2) or (is_amb1 and is_wt2):
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "WT"}
    if (is_wt1 and is_other_mut2) or (is_other_mut1 and is_wt2):
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "N/I"}

    if is_amb1 and is_amb2:
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "N/I"}
    if (is_amb1 and is_other_mut2) or (is_other_mut1 and is_amb2):
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "OTH"}

    if is_other_mut1 and is_other_mut2:
        return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "OTH"}

    return {"Detail": combine_original_statuses(mstat_5p, mstat_3p), "Simple": "UNK"}
