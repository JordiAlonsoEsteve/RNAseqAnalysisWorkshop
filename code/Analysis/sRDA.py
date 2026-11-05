import os
import re
import subprocess
import warnings

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D
from motifcompass.featureConstruction.k_mer import generate_kmer_count_vector
from motifcompass.featureConstruction.pwm import pwm_logsumexp_matrix
from motifcompass.model.cross_validation import tune_sparse_rda
from motifcompass.motif_de_novo.motif_construction import (
    generate_meme_file,
    generate_PWMs_from_A,
)
from motifcompass.sequenceConstruction.dna_regions import gff_to_fasta, read_fasta

SAMPLE_RE = re.compile(
    r"^(?P<temperature>-?\d+(?:\.\d+)?)_"
    r"(?P<treatment>.+)_"
    r"(?P<time>-?\d+(?:\.\d+)?)_"
    r"(?P<replicate>\d+)$"
)

def parse_sample_columns(cols):
    """Parse '<Temperature>_<Treatment>_<Time>_<Replicate>' into a DataFrame."""
    parts = pd.Series(pd.Index(cols).astype(str)).str.extract(SAMPLE_RE)
    bad = parts.isna().any(axis=1).to_numpy()
    if bad.any():
        raise ValueError(
            "Column names not in <Temperature>_<Treatment>_<Time>_<Replicate> form: "
            f"{list(np.asarray(cols)[bad])[:10]}"
        )
    parts["temperature"] = pd.to_numeric(parts["temperature"])
    parts["time"] = pd.to_numeric(parts["time"])
    parts["replicate"] = pd.to_numeric(parts["replicate"])
    return parts

def biplot(Y_centered, C, fitted_A, fitted_B,
           n_kmers=10, arrow_fill=0.65, size_range=(35, 250),
           markers=("o", "s", "^", "D", "v", "P")):
    """Biplot of the fitted RDA model (samples = points, k-mers/motifs = arrows)."""

    # --- order components by variance explained, highest first
    pct_comp = variance_explained(Y_centered, C, fitted_A, fitted_B)
    order = np.argsort(pct_comp)[::-1]
    A = fitted_A.values[:, order]
    B = fitted_B.values[order, :]

    # --- drop features with all-zero weights
    keep = np.any(A != 0, axis=1)
    A = A[keep]
    kmer_names = np.asarray(fitted_A.index)[keep]

    # --- sample metadata
    meta = parse_sample_columns(fitted_B.columns)
    treatment = meta["treatment"].to_numpy()
    temperature = meta["temperature"].to_numpy()
    time = meta["time"].to_numpy()

    # --- colours (treatment) and marker shapes (temperature)
    levels = list(pd.unique(treatment))
    cmap = plt.get_cmap("Set2" if len(levels) <= 8 else "tab20")
    colors = {h: cmap(i % cmap.N) for i, h in enumerate(levels)}

    temp_levels = sorted(pd.unique(temperature))
    marker_of = {t: markers[i % len(markers)] for i, t in enumerate(temp_levels)}

    # --- point sizes on a signed log scale so negative / zero times are allowed
    def tf(t):
        t = np.asarray(t, dtype=float)
        return np.sign(t) * np.log1p(np.abs(t))

    lt = tf(time)
    span = lt.max() - lt.min()

    def size_of(t):
        if span == 0:
            return np.full(np.shape(t), np.mean(size_range))
        return size_range[0] + (size_range[1] - size_range[0]) * (tf(t) - lt.min()) / span

    # --- coordinates
    S = B[[0, 1], :].T
    A2 = A[:, [0, 1]]
    length = np.linalg.norm(A2, axis=1)
    top = np.argsort(-length)[:n_kmers]
    arrows = A2[top] * (arrow_fill * np.abs(S).max() / length[top].max())

    # --- figure
    fig, ax = plt.subplots(figsize=(9, 6.5), dpi=150)
    ax.grid(color="#f0f0f0", linestyle="--", linewidth=1, zorder=0)
    ax.axhline(0, color="#d0d0d0", lw=1.2, zorder=1)
    ax.axvline(0, color="#d0d0d0", lw=1.2, zorder=1)

    for temp in temp_levels:
        m = temperature == temp
        ax.scatter(S[m, 0], S[m, 1], s=size_of(time[m]),
                   c=[colors[h] for h in treatment[m]],
                   marker=marker_of[temp],
                   alpha=0.85, edgecolor="white", lw=0.8, zorder=2)

    # --- arrows and labels
    texts = []
    arrow_color = "#A93226"
    halo = pe.withStroke(linewidth=2.5, foreground="white", alpha=0.9)

    for (x, y), j in zip(arrows, top):
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=arrow_color, lw=1.5,
                                    shrinkA=0, shrinkB=0, alpha=0.9),
                    zorder=3)
        txt = ax.text(x, y, kmer_names[j], color=arrow_color,
                      fontsize=10, fontweight="500", zorder=4)
        txt.set_path_effects([halo])
        texts.append(txt)
    ax.set_aspect("equal", adjustable="datalim")
    adjust_text(
        texts,
        x=S[:, 0], y=S[:, 1],
        target_x=arrows[:, 0], target_y=arrows[:, 1],
        ax=ax,
        force_text=(0.8, 1),
        force_static=(0.2, 0.3),
        force_pull=(0.01, 0.02),
        expand=(1.2, 1.35),
        max_move=(10, 10),
        ensure_inside_axes=True,
        expand_axes=True,
        iter_lim=1500,
        min_arrow_len=3,
        arrowprops=dict(arrowstyle="-", color=arrow_color, lw=0.7, alpha=0.65),
    )

    # --- legends: treatment (colour), temperature (shape), time (size)
    treatment_handles = [
        Line2D([], [], marker="o", ls="", color=colors[h], markersize=8,
               markeredgecolor="white", markeredgewidth=0.5, label=h)
        for h in levels
    ]
    leg1 = ax.legend(handles=treatment_handles, title="Treatment", frameon=False,
                     loc="upper left", title_fontproperties={"weight": "bold"})
    ax.add_artist(leg1)

    temp_handles = [
        Line2D([], [], marker=marker_of[t], ls="", color="#999999", markersize=8,
               markeredgecolor="white", markeredgewidth=0.5, label=f"{t:g} °C")
        for t in temp_levels
    ]
    leg2 = ax.legend(handles=temp_handles, title="Temperature", frameon=False,
                     loc="upper right", title_fontproperties={"weight": "bold"})
    ax.add_artist(leg2)

    ut = np.unique(time)
    show = np.unique(ut[np.linspace(0, len(ut) - 1, min(4, len(ut))).round().astype(int)])
    time_handles = [
        Line2D([], [], marker="o", ls="", color="#999999",
               markersize=np.sqrt(float(size_of(np.array([t]))[0])),
               markeredgecolor="white", markeredgewidth=0.5, label=f"{t:g} h")
        for t in show
    ]
    ax.legend(handles=time_handles, title="Time", frameon=False,
              loc="lower left", labelspacing=1.2,
              title_fontproperties={"weight": "bold"})

    # --- axes
    total = pct_comp.sum()
    for axis_idx, setter in enumerate((ax.set_xlabel, ax.set_ylabel)):
        c = order[axis_idx]
        setter(f"Component {c + 1} (R² = {pct_comp[c]:.1f}%, "
               f"Prop. Var = {100 * pct_comp[c] / total:.2f}%)",
               fontsize=11, fontweight="bold", color="#444444")

    ax.tick_params(axis="both", colors="#666666", labelsize=9, length=4)
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_color("#bbbbbb")
        ax.spines[s].set_linewidth(1.2)

    fig.tight_layout()
    return ax

def generate_promoter(gff_file, fasta_file, output_fasta, upstream=1000, downstream=100):
    promoter_dict = gff_to_fasta(fasta_file, gff_file, output_fasta, upstream=upstream, downstream=downstream)
    return promoter_dict

def generate_X_pwm(promoter_dict, meme_file, pseudocount=0.01, both_strands=True):
    X_pwm = pwm_logsumexp_matrix(promoter_dict, meme_file, pseudocount=pseudocount, both_strands=both_strands)
    return X_pwm

def variance_explained(Y_centered, C, fitted_A, fitted_B):
    """Compute variance explained by each component."""
    Z = C.values @ fitted_A.values
    tss = (Y_centered.values ** 2).sum()
    pct_comp = []
    for i in range(fitted_A.shape[1]):
        fitted_values = np.outer(Z[:, i], fitted_B.values[i, :])
        rss = ((Y_centered.values - fitted_values) ** 2).sum()
        pct_comp.append(100 * (tss - rss) / tss)
    return np.asarray(pct_comp)

def clean_gene_ids(index):
    """Drop '::coords'"""
    ids = pd.Index(index).astype(str)
    ids = ids.str.split("::").str[0]
    return ids

def load_Y(path):
    """Genes x samples table (columns <Temp>_<Treatment>_<Time>_<Rep>)"""
    Y = pd.read_csv(path, index_col="gene_id")
    parse_sample_columns(Y.columns)  # fail early if a column name is malformed
    
    return Y

def main(Y_data_url, X_type, out_dir, study_name):
    
    os.makedirs(out_dir, exist_ok=True)

    # --- promoters
    promoter_fasta = f"{out_dir}/promoters.fa"
    generate_promoter(
        gff_file="data/genomic_data/TAIR_12/annotation.gff",
        fasta_file="data/genomic_data/TAIR_12/col-0_keygene.fasta",
        output_fasta=promoter_fasta,
        upstream=500,
        downstream=200,
    )
    promoters = read_fasta(promoter_fasta)

    # --- X matrix (genes x features)
    meme_file = "data/binding/jaspar_athaliana.meme"

    if X_type == "kmer":
        X = generate_kmer_count_vector(promoters=promoters, k=6, both_strands=True)
    elif X_type == "pwm":
        X = generate_X_pwm(promoters, meme_file, pseudocount=0.01, both_strands=True)
    else:
        raise ValueError(f"Unknown X_type: {X_type}")

    X.index = clean_gene_ids(X.index)
    X = X[~X.index.duplicated(keep="first")]  # one row per gene

    # --- Y matrix (genes x samples)
    Y = load_Y(Y_data_url)
    Y.index = clean_gene_ids(Y.index)

    shared = X.index.intersection(Y.index)
    if len(shared) < 100:
        raise ValueError(
            f"Only {len(shared)} genes shared between X and Y. "
            "Check that the gene ID formats of promoters.fa and the counts table match."
        )
    X, Y = X.loc[shared], Y.loc[shared]
    print(f"X shape: {X.shape}, Y shape: {Y.shape}")

    # --- tune and fit sparse RDA
    best = tune_sparse_rda(
        X.values, Y.values,
        storage_dir=out_dir, study_name=f"{study_name}_{X_type}_rda",
        alpha=(0.01, 0.1), n_trials=20,
        n_components=(2, 10), n_splits=3,
        standardize_x=(X_type == "kmer"),
        verbose=False, max_iter=100, tol=1e-4,
    )

    comp_names = [f"Comp_{i + 1}" for i in range(best.A.shape[1])]
    A_df = pd.DataFrame(best.A, index=X.columns, columns=comp_names)
    B_df = pd.DataFrame(best.B, index=comp_names, columns=Y.columns)
    Y_centered = Y - Y.mean(axis=0)

    A_df.to_csv(f"{out_dir}/A_loadings.tsv", sep="\t")
    B_df.to_csv(f"{out_dir}/B_sample_scores.tsv", sep="\t")

    # --- k-mer specific: standardise X, build motifs, compare to JASPAR
    if X_type == "kmer":
        X = (X - X.mean(axis=0)) / X.std(axis=0)
        pwms = generate_PWMs_from_A(A_df, ic_threshold=10, signed=True, threshold_increase=0.75)
        motifs_meme = f"{out_dir}/motifs.meme"
        generate_meme_file(pwms, motifs_meme)
        subprocess.run(
            ["tomtom", "-oc", f"{out_dir}/tomtom", "-thresh", "0.05",
             motifs_meme, meme_file, "-verbosity", "2"],
            check=True,
        )

    # --- biplot
    ax = biplot(Y_centered, X, A_df, B_df, n_kmers=20)
    ax.figure.savefig(f"{out_dir}/biplot_{X_type}.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main(X_type="pwm", Y_data_url="data/hormone_temp/processed/log2FC.csv", out_dir="data/hormone_temp/analysis_rda", study_name="hormone")

