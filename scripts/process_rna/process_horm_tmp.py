import re
from pathlib import Path
import pandas as pd

# -------------------- EDIT THESE SETTINGS --------------------
# Relative paths are resolved from your current working directory.
TSV_FILE = Path("data/hormone_temp/salmon.merged.gene_counts.tsv")
EXCEL_FILE = Path("data/hormone_temp/Metadata_bigexperiment_hormones.xlsx")
EXCEL_SHEET = 0  # First worksheet; alternatively use a name, e.g. "Sheet1"
OUTPUT_FILE = Path("data/hormone_temp/counts_genes_x_samples.tsv")

# WARNING: samples whose Treatment CONTAINS any of these words (case-insensitive)
# are DROPPED from the output. This is substring matching, so a treatment such as
# "room_temp_control" or "pretreated_room" would also be removed. Check the printed
# list of dropped samples below to confirm only the intended ones are excluded.
DROP_TREATMENT_KEYWORDS = ["untreated", "room"]
# -------------------------------------------------------------
GENE_PREFIX = re.compile(r"^gene-TAIR12_TAIR12_")

METADATA_COLUMNS = [
    "BMKgene_id",
    "lab_internal_sample_id",
    "Temperature",
    "Treatment",
    "Timepoint",
    "Replicate",
]


def clean_text(series):
    """Trim whitespace while preserving missing values."""
    return series.astype("string").str.strip().replace("", pd.NA)


def check_unique_ids(series, label):
    if series.isna().any():
        raise ValueError(f"{label} contains missing or blank identifiers.")
    duplicates = series[series.duplicated(keep=False)].unique().tolist()
    if duplicates:
        raise ValueError(f"Duplicate {label} values: {duplicates[:10]}")


def main():
    genes = pd.read_csv(TSV_FILE, sep="\t", dtype={"gene_id": "string"})
    metadata = pd.read_excel(EXCEL_FILE, sheet_name=EXCEL_SHEET)

    # Clean headers and reject ambiguous column names.
    for table, label in [(genes, "TSV"), (metadata, "Excel")]:
        table.columns = table.columns.astype(str).str.strip()
        if table.columns.duplicated().any():
            raise ValueError(f"{label} contains duplicate column names.")

    genes["gene_id"] = genes["gene_id"].astype("string").str.strip()
    genes["gene_id"] = genes["gene_id"].str.replace(GENE_PREFIX, "", regex=True)
    genes["gene_name"] = genes["gene_name"].astype("string").str.strip()
    genes["gene_name"] = genes["gene_name"].str.replace(GENE_PREFIX, "", regex=True)

    if genes["gene_id"].isna().any() or genes["gene_id"].duplicated().any():
        dupes = genes.loc[genes["gene_id"].duplicated(keep=False), "gene_id"].unique().tolist()
        raise ValueError(f"gene_id has missing or duplicate values after stripping: {dupes[:10]}")

        
    if "Timepoint" not in metadata.columns:
        metadata = metadata.rename(columns={"Timepoint_in_hours": "Timepoint"})

    missing_gene_columns = {"gene_id", "gene_name"} - set(genes.columns)
    if missing_gene_columns:
        raise ValueError(f"Missing TSV columns: {sorted(missing_gene_columns)}")

    missing_metadata = set(METADATA_COLUMNS) - set(metadata.columns)
    if missing_metadata:
        raise ValueError(f"Missing Excel columns: {sorted(missing_metadata)}")

    metadata = metadata[METADATA_COLUMNS].copy()
    for column in ["BMKgene_id", "lab_internal_sample_id", "Treatment"]:
        metadata[column] = clean_text(metadata[column])
    genes["gene_id"] = clean_text(genes["gene_id"])

    check_unique_ids(genes["gene_id"], "gene_id")
    check_unique_ids(metadata["BMKgene_id"], "BMKgene_id")

    # Convert numeric metadata, including decimal commas such as '0,5'.
    for column in ["Temperature", "Timepoint", "Replicate"]:
        values = clean_text(metadata[column]).str.replace(",", ".", regex=False)
        metadata[column] = pd.to_numeric(values, errors="raise")

    # Every TSV column except these two annotation columns is a sample.
    sample_columns = [
        column for column in genes.columns
        if column not in ["gene_id", "gene_name"]
    ]
    if not sample_columns:
        raise ValueError("No sample measurement columns were found in the TSV.")

    missing_samples = sorted(set(sample_columns) - set(metadata["BMKgene_id"]))
    if missing_samples:
        raise ValueError(
            "TSV sample IDs missing from Excel metadata: "
            + ", ".join(missing_samples[:20])
            + ". Every TSV sample must have a matching BMKgene_id."
        )

    extra_samples = sorted(set(metadata["BMKgene_id"]) - set(sample_columns))
    if extra_samples:
        print(f"Note: ignoring {len(extra_samples)} Excel samples absent from the TSV.")

    # Align metadata to TSV sample order.
    aligned_metadata = metadata.set_index("BMKgene_id").loc[sample_columns]

    # WARNING: drop samples whose Treatment contains any keyword in
    # DROP_TREATMENT_KEYWORDS (case-insensitive substring match).
    #pattern = "|".join(re.escape(word) for word in DROP_TREATMENT_KEYWORDS)
    #drop_mask = aligned_metadata["Treatment"].str.contains(pattern, case=False, na=False)
    #dropped = aligned_metadata.loc[drop_mask, "Treatment"]
    #if not dropped.empty:
    #    print(
    #        f"WARNING: dropping {len(dropped)} sample(s) whose Treatment contains "
    #        f"{DROP_TREATMENT_KEYWORDS}:"
    #    )
    #    print(dropped.to_string())
    #sample_columns = [c for c in sample_columns if c not in dropped.index]
    if not sample_columns:
        raise ValueError("All samples were dropped by the Treatment filter.")
    aligned_metadata = aligned_metadata.loc[sample_columns]

    # Build new column names (e.g. 21_PAMP_0.5_1) in TSV sample order.
    new_names = (
        aligned_metadata["Temperature"].map("{:g}".format) + "_"
        + aligned_metadata["Treatment"].astype(str) + "_"
        + aligned_metadata["Timepoint"].map("{:g}".format) + "_"
        + aligned_metadata["Replicate"].map("{:g}".format)
    )
    if new_names.duplicated().any():
        raise ValueError(
            f"New column names are not unique: {new_names[new_names.duplicated(keep=False)].unique().tolist()[:10]}"
        )

    # Genes stay as rows, samples as columns (no transposing).
    out = genes[["gene_id", "gene_name"] + sample_columns].copy()
    out[sample_columns] = out[sample_columns].apply(pd.to_numeric, errors="raise")
    out = out.rename(columns=dict(zip(sample_columns, new_names)))

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print(f"Saved: {OUTPUT_FILE.resolve()}")
    print(f"Genes: {len(out)} | Samples kept: {len(sample_columns)}")
    print(out.iloc[:5, :6].to_string(index=False))


if __name__ == "__main__":
    main()