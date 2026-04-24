"""Fragment size distribution plot.

Mirrors plot_size_distribution() in R/plot_size_distribution.R.
Requires matplotlib; seaborn is used for density curves if available,
otherwise scipy.stats.gaussian_kde is used as a fallback.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import pandas as pd


# Nucleosome peak positions (bp) used as reference lines
_NUC_PEAKS = {
    "mono":  167,
    "di":    334,
    "tri":   500,
}


def plot_size_distribution(
    df_fragments: pd.DataFrame,
    size_col: str = "Fragment_Size",
    col_z: Optional[str] = "Fragment_Status_Simple",
    vals_z: Optional[Sequence[str]] = None,
    show_histogram: bool = False,
    show_density: bool = True,
    x_limits: Tuple[float, float] = (0, 600),
    histogram_binwidth: float = 5,
    histo_args: Optional[Dict[str, Any]] = None,
    density_args: Optional[Dict[str, Any]] = None,
    colors_z: Optional[Sequence[str]] = None,
    show_nuc_peaks: bool = True,
    title: Optional[str] = None,
    output_path: Optional[str] = None,
) -> plt.Figure:
    """Plot the fragment size distribution, optionally grouped and/or saved.

    Parameters
    ----------
    df_fragments:
        DataFrame produced by ``run_fragmentomYcs``.
    size_col:
        Column name containing fragment lengths (numeric).
    col_z:
        Column name used for grouping. ``None`` for no grouping.
    vals_z:
        Subset of ``col_z`` values to display. ``None`` = all groups.
    show_histogram:
        Add histogram bars.
    show_density:
        Add kernel-density curves.
    x_limits:
        (min, max) range for the x-axis.
    histogram_binwidth:
        Bin width in base pairs.
    histo_args:
        Extra kwargs forwarded to ``ax.hist()``.
    density_args:
        Extra kwargs forwarded to the KDE plot call.
    colors_z:
        Sequence of hex/named colours (one per group, in order). ``None``
        uses the matplotlib default colour cycle.
    show_nuc_peaks:
        Draw dashed vertical lines at the mono/di/tri-nucleosome peak positions.
    title:
        Plot title.  ``None`` generates a default title.
    output_path:
        If given, the figure is saved to this path and the function returns
        ``None``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if histo_args is None:
        histo_args = {}
    if density_args is None:
        density_args = {}

    if not show_histogram and not show_density:
        raise ValueError("At least one of show_histogram or show_density must be True.")

    if size_col not in df_fragments.columns:
        raise ValueError(f"Column '{size_col}' not found in the dataframe.")
    if not pd.api.types.is_numeric_dtype(df_fragments[size_col]):
        raise ValueError(f"Column '{size_col}' must be numeric.")
    if col_z is not None and col_z not in df_fragments.columns:
        raise ValueError(f"Column '{col_z}' not found in the dataframe.")
    if col_z is None and vals_z is not None:
        raise ValueError("If col_z is None, vals_z must also be None.")

    df = df_fragments.copy()

    # ---- Filter to requested groups
    if col_z is not None and vals_z is not None:
        df = df[df[col_z].isin(vals_z)]

    # ---- Build group list
    if col_z is None:
        groups = [("All", df)]
    else:
        order = list(vals_z) if vals_z is not None else sorted(df[col_z].dropna().unique())
        groups = [(g, df[df[col_z] == g]) for g in order]

    # ---- Assign colours
    if colors_z is not None:
        if len(colors_z) < len(groups):
            raise ValueError(
                f"Provided {len(colors_z)} colour(s) but need {len(groups)} for: "
                + ", ".join(g for g, _ in groups)
            )
        colour_map = {g: colors_z[i] for i, (g, _) in enumerate(groups)}
    else:
        prop_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        colour_map = {g: prop_cycle[i % len(prop_cycle)] for i, (g, _) in enumerate(groups)}

    fig, ax = plt.subplots()

    for group_name, df_group in groups:
        sizes = df_group[size_col].dropna()
        if show_density and len(sizes) < 2:
            warnings.warn(
                f"Group '{group_name}' has fewer than 2 observations — density skipped.",
                stacklevel=2,
            )
            show_density_this = False
        else:
            show_density_this = show_density

        n = len(sizes)
        label = f"{group_name} (N={n})"
        colour = colour_map[group_name]

        # ---- Histogram
        if show_histogram:
            bins = _make_bins(x_limits, histogram_binwidth)
            h_defaults = {"alpha": 0.5, "density": show_density_this}
            h_kw = {**h_defaults, **histo_args, "color": colour, "label": label if not show_density_this else "_nolegend_"}
            ax.hist(sizes, bins=bins, **h_kw)

        # ---- Density
        if show_density_this:
            d_defaults = {"linewidth": 1.5}
            d_kw = {**d_defaults, **density_args, "color": colour, "label": label}
            _plot_kde(ax, sizes, x_limits, **d_kw)

    # ---- Nucleosome peaks
    if show_nuc_peaks:
        for peak_name, peak_pos in _NUC_PEAKS.items():
            if x_limits[0] <= peak_pos <= x_limits[1]:
                ax.axvline(peak_pos, color="grey", linestyle="--", linewidth=0.8,
                           label=f"{peak_name} ({peak_pos} bp)")

    ax.set_xlim(x_limits)
    ax.set_xlabel("Fragment size (bp)")
    ax.set_ylabel("Density" if show_density else "Count")
    ax.set_title(title or "Fragment size distribution")
    ax.legend()

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        return None

    return fig


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _make_bins(x_limits: Tuple[float, float], binwidth: float) -> List[float]:
    """Return bin edges covering *x_limits* with step *binwidth*."""
    import numpy as np
    return list(np.arange(x_limits[0], x_limits[1] + binwidth, binwidth))


def _plot_kde(ax: plt.Axes, data: pd.Series, x_limits: Tuple[float, float], **kwargs) -> None:
    """Draw a KDE curve on *ax*.  Uses scipy if available, else numpy."""
    import numpy as np
    x = np.linspace(x_limits[0], x_limits[1], 512)
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(data)
        y = kde(x)
    except ImportError:
        # numpy-based fallback (less smooth)
        counts, edges = np.histogram(data, bins=100, range=x_limits, density=True)
        centres = (edges[:-1] + edges[1:]) / 2
        y = np.interp(x, centres, counts)
    ax.plot(x, y, **kwargs)
