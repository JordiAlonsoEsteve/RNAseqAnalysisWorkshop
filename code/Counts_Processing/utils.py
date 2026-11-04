import pandas as pd
import numpy as np
import warnings

 
SAMPLE_FIELDS = ["temperature", "treatment", "time", "replicate"]
 
def parse_sample_names(columns, sep="_") -> pd.DataFrame:
    """Split '<Temperature>_<Treatment>_<Time>_<Replicate>' column names into a metadata table."""
    columns = pd.Index(columns)
    if not columns.is_unique:
        raise ValueError("Column names must be unique.")
    parts = [str(c).split(sep) for c in columns]
    bad = [c for c, p in zip(columns, parts) if len(p) != len(SAMPLE_FIELDS)]
    if bad:
        raise ValueError(
            f"Columns must have the form <Temperature>{sep}<Treatment>{sep}<Time>{sep}<Replicate>. "
            f"Offending columns: {bad[:5]}"
        )
    return pd.DataFrame(parts, index=columns, columns=SAMPLE_FIELDS)

def median_of_ratios_normalization(counts: pd.DataFrame) -> pd.DataFrame:

    # Step 1: Compute geometric mean for each gene
    def geo_mean_nonzero(x):
        x = x[x > 0]
        if len(x) == 0:
            return np.nan
        return np.exp(np.mean(np.log(x)))

    geometric_means = counts.apply(geo_mean_nonzero, axis=1)

    if geometric_means.isna().any():
        n_drop = geometric_means.isna().sum()
        warnings.warn(f"Dropping {n_drop} genes with zero counts in all samples")
        keep = geometric_means.notna()
        counts = counts.loc[keep]
        geometric_means = geometric_means.loc[keep]

    # Step 2: Compute ratios of counts to geometric means
    # ratios = counts.div(geometric_means, axis=0)
    ratios = counts.div(geometric_means.replace(0, np.nan), axis=0)

    # Step 3: Compute median ratio per sample: size factors
    size_factors = ratios.median(axis=0, skipna=True)
    zero_sf = size_factors[size_factors == 0]
    if not zero_sf.empty:
        warnings.warn(
            f"Size factor is zero for samples: {list(zero_sf.index)}, dropping"
        )
        counts = counts.drop(columns=zero_sf.index)
        ratios = ratios.drop(columns=zero_sf.index)
        size_factors = size_factors.drop(index=zero_sf.index)

    # Step 4: Normalize counts
    normalized_counts = counts.div(size_factors, axis=1)

    return normalized_counts

def process_counts_log2FC(
    counts: pd.DataFrame,
    control_temperature: str = "21C",
    control_treatment: str = "control",
    average_treatments: bool = False,
    clipping: float = 4,
    pseudocount: float = 1,
    sep: str = "_",
) -> pd.DataFrame:
    """
    log2 fold change of every sample w.r.t. a common reference.
 
    Columns must be '<Temperature>_<Treatment>_<Time>_<Replicate>'.
    The reference samples are those with Temperature == `control_temperature`
    AND Treatment == `control_treatment`. The same reference condition is used for
    all other temperatures and treatments, taken at the same time point as the
    sample being compared (reference samples are averaged over replicates only).
 
    log2FC = clip( log2((K + pc) / (ref_t + pc)), -clipping, clipping )
 
    Pass (median-of-ratios) normalized counts. The reference samples themselves are
    not returned. If `average_treatments` is True, replicates of the remaining
    samples are averaged first and columns are named '<Temperature>_<Treatment>_<Time>'.
    Returns genes x non-reference samples.
    """
    counts = counts.dropna(axis=0, how="any")
    meta = parse_sample_names(counts.columns, sep)
 
    control_temperature = str(control_temperature)
    is_ref = (meta["temperature"] == control_temperature) & (
        meta["treatment"] == control_treatment
    )
    if not is_ref.any():
        raise ValueError(
            f"No reference samples found (temperature == '{control_temperature}' "
            f"and treatment == '{control_treatment}')."
        )
    if is_ref.all():
        raise ValueError("Only reference samples found; nothing to compare.")
 
    # Reference: mean over replicates, one value per gene and time point
    reference = (
        counts.loc[:, is_ref].astype(float).T.groupby(meta.loc[is_ref, "time"]).mean().T
    )
 
    # Remaining samples (optionally averaged over replicates)
    other_meta = meta[~is_ref]
    other_counts = counts.loc[:, ~is_ref].astype(float)
    if average_treatments:
        group = other_meta[["temperature", "treatment", "time"]].agg(sep.join, axis=1)
        other_counts = other_counts.T.groupby(group).mean().T
        other_meta = (
            other_meta.assign(group=group)
            .drop_duplicates("group")
            .set_index("group")
            .loc[other_counts.columns]
        )
 
    # Match each sample to the reference at its own time point
    missing = sorted(set(other_meta["time"]) - set(reference.columns))
    if missing:
        raise ValueError(f"No reference samples at time(s): {missing}")
    ref_aligned = reference[other_meta["time"].values]
    ref_aligned.columns = other_counts.columns
 
    lfc = np.log2((other_counts + pseudocount) / (ref_aligned + pseudocount))
    lfc = lfc.clip(-clipping, clipping)
 
    lfc = lfc.sort_index()
    lfc = lfc.loc[:, lfc.std() != 0]  # drop constant columns
    lfc = lfc.dropna()
 
    print(f"Final shape, {lfc.shape}")
    return lfc

if __name__ == "__main__":
    raw_pdamp = pd.read_csv("data/PDAMP/counts_genes_x_samples.tsv", index_col=0, sep="\t")
    # drop gene_name column if present
    if "gene_name" in raw_pdamp.columns:
        raw_pdamp = raw_pdamp.drop(columns="gene_name")
    # drop columns that have untreated as treatment
    raw_pdamp = raw_pdamp.loc[:, ~raw_pdamp.columns.str.contains("_untreated_")]
    
    normalized_pdamp_data = median_of_ratios_normalization(raw_pdamp)

    processed_pdamp_data = process_counts_log2FC(normalized_pdamp_data, control_temperature = "21", control_treatment = "mock")
    # save 
    processed_pdamp_data.to_csv("data/PDAMP/processed/log2FC.csv")