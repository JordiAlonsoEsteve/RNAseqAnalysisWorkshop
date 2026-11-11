import re
import sys
sys.path.append("code")  # for utils.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from Counts_Processing.utils import parse_sample_names
# First 8 colors: reference palette (fixed order). Beyond 8, tab20 / turbo are used.
BASE_COLORS = [
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
    "#e87ba4", "#008300", "#6250d6", "#e34948",
]
# Marker shapes in order of the levels of `marker_by` ("o" and "x" for 2 temperatures)
MARKERS = ["o", "x", "s", "^", "D", "v"]
UNFILLED = {"x", "+", "1", "2", "3", "4"}  # drawn with lines only: color goes on the line

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"


def compute_pca(lfc: pd.DataFrame, standardize: bool = False, n_components: int = 10):
    """
    PCA of samples, with genes as features.

    lfc: genes x samples (e.g. output of process_counts_log2FC).
    The mean of each gene (across samples) is ALWAYS removed.
    If standardize=True, each gene is also divided by its standard deviation,
    so every gene contributes equally (otherwise high-variance genes dominate).

    Returns (scores, explained_variance_ratio, loadings):
      scores    : samples x PCs DataFrame
      explained : Series with the fraction of variance per PC
      loadings  : genes x PCs DataFrame
    """
    X = lfc.T.astype(float)  # samples x genes

    # Center each gene (always); divide by its std only if standardize=True
    X_scaled = StandardScaler(with_mean=True, with_std=standardize).fit_transform(X)

    k = min(n_components, X.shape[0], X.shape[1])
    pca = PCA(n_components=k)
    scores_arr = pca.fit_transform(X_scaled)  # samples x PCs

    pcs = [f"PC{i + 1}" for i in range(k)]
    scores = pd.DataFrame(scores_arr, index=X.index, columns=pcs)
    explained = pd.Series(pca.explained_variance_ratio_, index=pcs)
    loadings = pd.DataFrame(pca.components_.T, index=X.columns, columns=pcs)
    return scores, explained, loadings


def _levels(series: pd.Series) -> list:
    """Unique levels, sorted numerically if they all look like numbers, else in order of appearance."""
    levels = list(dict.fromkeys(series))
    try:
        return sorted(levels, key=lambda v: float(v))
    except ValueError:
        return levels


def _numeric(series: pd.Series) -> np.ndarray:
    """Numeric value of labels like '1h', '24h', '0.5' (the first number in the string).
    Use a common unit for all time points. Falls back to the rank of the level."""
    nums = series.astype(str).str.extract(r"(\d+(?:\.\d+)?)")[0].astype(float)
    if nums.isna().any():
        order = {lvl: i + 1 for i, lvl in enumerate(_levels(series))}
        return series.map(order).astype(float).to_numpy()
    return nums.to_numpy()


def _palette(n: int) -> list:
    if n <= len(BASE_COLORS):
        return BASE_COLORS[:n]
    if n <= 20:
        return [plt.get_cmap("tab20")(i) for i in range(n)]
    return [plt.get_cmap("turbo")(v) for v in np.linspace(0.05, 0.95, n)]


def plot_pca(
    lfc: pd.DataFrame,
    standardize: bool = False,
    pcs=(1, 2),
    color_by: str = "treatment",
    marker_by: str = "temperature",
    size_by: str | None = "time",
    size_range=(40, 220),
    colors: dict | None = None,
    annotate: str | None = None,
    label_groups: bool = False,
    title: str | None = None,
    ax=None,
    sep: str = "_",
):
    """
    PCA scatter of samples (columns named '<Temperature>_<Treatment>_<Time>_<Replicate>').

      color  -> `color_by`  (any number of levels; pass `colors={level: color}` to fix them)
      shape  -> `marker_by` ("o", "x", then "s", "^", ...)
      size   -> `size_by`   (numeric value of the label, log-scaled; None = constant size)

    `annotate` writes one of the fields next to every point; `label_groups` writes each
    color group's name next to its points.
    """
    scores, explained, _ = compute_pca(lfc, standardize, n_components=max(pcs))
    meta = parse_sample_names(scores.index, sep)

    color_levels = _levels(meta[color_by])
    marker_levels = _levels(meta[marker_by])
    if len(marker_levels) > len(MARKERS):
        raise ValueError(f"{marker_by} has too many levels ({len(marker_levels)}) for distinct markers.")
    color_of = dict(zip(color_levels, _palette(len(color_levels))))
    if colors:
        color_of.update(colors)
    marker_of = dict(zip(marker_levels, MARKERS))

    # Marker size from the numeric value of `size_by` (log scale when all values > 0)
    if size_by is not None:
        vals = _numeric(meta[size_by])
        v = np.log(vals) if (vals > 0).all() else vals
        span = v.max() - v.min()
        frac = (v - v.min()) / span if span > 0 else np.full_like(v, 0.5)
        sizes = pd.Series(size_range[0] + (size_range[1] - size_range[0]) * frac, index=meta.index)
    else:
        sizes = pd.Series(np.mean(size_range), index=meta.index)

    created = ax is None
    if created:
        fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=150)
    fig = ax.figure
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    px, py = f"PC{pcs[0]}", f"PC{pcs[1]}"
    for ml in marker_levels:
        idx = meta.index[meta[marker_by] == ml]
        mk = marker_of[ml]
        c = [color_of[lvl] for lvl in meta.loc[idx, color_by]]
        if mk in UNFILLED:  # "x": drawn as lines, so color the lines
            kw = dict(edgecolors="face", linewidths=1.8)
        else:  # filled: thin surface ring keeps overlapping points readable
            kw = dict(edgecolors=SURFACE, linewidths=0.5)
        ax.scatter(
            scores.loc[idx, px], scores.loc[idx, py],
            s=sizes[idx], c=c, marker=mk, zorder=3, **kw, alpha=0.9
        )

    if label_groups:
        halo = pe.withStroke(linewidth=2.5, foreground=SURFACE)
        for cl in color_levels:
            idx = meta.index[meta[color_by] == cl]
            xy = scores.loc[idx, [px, py]]
            nearest = ((xy - xy.mean()) ** 2).sum(axis=1).idxmin()  # member closest to centroid
            ax.annotate(
                cl, (xy.loc[nearest, px], xy.loc[nearest, py]),
                xytext=(0, 11), textcoords="offset points",
                color=INK, fontsize=8, fontweight="bold", ha="center", va="bottom",
                zorder=5, path_effects=[halo],
            )

    if annotate is not None:
        for name in scores.index:
            ax.annotate(
                meta.loc[name, annotate], (scores.loc[name, px], scores.loc[name, py]),
                xytext=(6, 5), textcoords="offset points", fontsize=7, color=INK_2, zorder=4,
            )

    # Recessive chrome: hairline grid, muted axes, zero lines
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    ax.axhline(0, color=AXIS, linewidth=0.8, zorder=1)
    ax.axvline(0, color=AXIS, linewidth=0.8, zorder=1)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=8, length=3)

    ax.set_xlabel(f"{px} ({explained[px] * 100:.1f}% variance)", color=INK_2, fontsize=9)
    ax.set_ylabel(f"{py} ({explained[py] * 100:.1f}% variance)", color=INK_2, fontsize=9)
    mode = "standardized" if standardize else "centered"
    ax.set_title(
        title or f"PCA of log$_2$FC samples ({mode} genes)",
        loc="left", color=INK, fontsize=11, fontweight="bold",
    )

    # One legend with three sections: color, shape, size (text stays in ink)
    entries = [(Line2D([], [], linestyle=""), color_by, True)]
    entries += [
        (Line2D([], [], marker="o", linestyle="", markersize=8, color=color_of[c],
                markeredgecolor=SURFACE), str(c), False)
        for c in color_levels
    ]
    if len(marker_levels) > 1:
        entries.append((Line2D([], [], linestyle=""), marker_by, True))
        for m in marker_levels:
            mk = marker_of[m]
            entries.append((
                Line2D([], [], marker=mk, linestyle="", markersize=8, color=INK_2,
                       markeredgecolor=INK_2 if mk in UNFILLED else SURFACE,
                       markeredgewidth=1.8 if mk in UNFILLED else 1.0),
                str(m), False,
            ))
    if size_by is not None:
        entries.append((Line2D([], [], linestyle=""), size_by, True))
        uniq = _levels(meta[size_by])
        show = uniq if len(uniq) <= 5 else [uniq[i] for i in np.unique(np.linspace(0, len(uniq) - 1, 4).round().astype(int))]
        for lvl in show:
            s = sizes[meta.index[meta[size_by] == lvl][0]]
            entries.append((
                Line2D([], [], marker="o", linestyle="", markersize=np.sqrt(s), color="#a9a8a2",
                       markeredgecolor=SURFACE),
                str(lvl), False,
            ))

    leg = ax.legend(
        handles=[e[0] for e in entries], labels=[e[1] for e in entries],
        loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8,
        labelspacing=0.8, handlelength=1.2, ncol=1,
    )
    for text, (_, _, is_header) in zip(leg.get_texts(), entries):
        text.set_color(INK_2 if is_header else INK)
        if is_header:
            text.set_fontweight("bold")

    if created:
        fig.tight_layout()
    return ax


if __name__ == "__main__":
    lfc = pd.read_csv("data/hormone_temp/processed/log2FC.csv", index_col=0)
    # drop if + in the column name
    lfc = lfc.loc[:, ~lfc.columns.str.contains(r"\+")]
    # drop columns that have 28C
    #lfc = lfc.loc[:, ~lfc.columns.str.contains("28C")]
    ax = plot_pca(lfc, standardize=False)
    import os
    os.makedirs("data/hormone_temp/analysis_pca", exist_ok=True)
    plt.savefig("data/hormone_temp/analysis_pca/pca_example_no_int.png", dpi=300)
    print("saved pca_example.png")

    # Now PC 3 and 4
    ax = plot_pca(lfc, standardize=False, pcs=(3, 4))
    plt.savefig("data/hormone_temp/analysis_pca/pca_example_no_int_PC3_4.png", dpi=300)

    # Now 5 and 6
    ax = plot_pca(lfc, standardize=False, pcs=(5, 6))
    plt.savefig("data/hormone_temp/analysis_pca/pca_example_no_int_PC5_6.png", dpi=300)