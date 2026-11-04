import re
from pathlib import Path
import pandas as pd


# -------------------- EDIT THESE SETTINGS --------------------
TSV_FILE = Path("data/PDAMP/salmon.merged.gene_counts_PAMPDAMP.tsv")
EXCEL_FILE = Path("data/PDAMP/Metadata_PAMPDAMP_SimpleVersion.xlsx")
EXCEL_SHEET = 0
OUTPUT_FILE = Path("data/PDAMP/counts_genes_x_samples.tsv")
KEEP_GENE_NAME = True  # keep gene_name as a second column next to gene_id
# -------------------------------------------------------------

GENE_PREFIX = re.compile(r"^gene-TAIR12_TAIR12_")
METADATA_COLUMNS = ["BMKgene_id", "Temperature", "Treatment", "Timepoint", "Replicate"]


def fmt(value):
    """21.0 -> '21', 0.5 -> '0.5', -24.0 -> '-24'."""
    return f"{value:g}"


def main():
    genes = pd.read_csv(TSV_FILE, sep="\t", dtype={"gene_id": "string"})
    metadata = pd.read_excel(EXCEL_FILE, sheet_name=EXCEL_SHEET)

    genes.columns = genes.columns.astype(str).str.strip()
    metadata.columns = metadata.columns.astype(str).str.strip()
    genes["gene_id"] = genes["gene_id"].astype("string").str.strip()
    genes["gene_id"] = genes["gene_id"].str.replace(GENE_PREFIX, "", regex=True)
    genes["gene_name"] = genes["gene_name"].astype("string").str.strip()
    genes["gene_name"] = genes["gene_name"].str.replace(GENE_PREFIX, "", regex=True)

    if genes["gene_id"].isna().any() or genes["gene_id"].duplicated().any():
        dupes = genes.loc[genes["gene_id"].duplicated(keep=False), "gene_id"].unique().tolist()
        raise ValueError(f"gene_id has missing or duplicate values after stripping: {dupes[:10]}")
    

    if "Timepoint" not in metadata.columns:
        metadata = metadata.rename(columns={"Timepoint_in_hours": "Timepoint"})

    missing = set(METADATA_COLUMNS) - set(metadata.columns)
    if missing:
        raise ValueError(f"Missing Excel columns: {sorted(missing)}")

    metadata = metadata[METADATA_COLUMNS].copy()
    metadata["BMKgene_id"] = metadata["BMKgene_id"].astype("string").str.strip()
    metadata["Treatment"] = metadata["Treatment"].astype("string").str.strip()

    # Numeric metadata, tolerating decimal commas such as '0,5'
    for col in ["Temperature", "Timepoint", "Replicate"]:
        values = metadata[col].astype("string").str.strip().str.replace(",", ".", regex=False)
        metadata[col] = pd.to_numeric(values, errors="raise")

    # Build the new column name per sample
    metadata["new_name"] = (
        metadata["Temperature"].map(fmt) + "_"
        + metadata["Treatment"] + "_"
        + metadata["Timepoint"].map(fmt) + "_"
        + metadata["Replicate"].map(fmt)
    )
    name_map = metadata.set_index("BMKgene_id")["new_name"]

    # Checks
    sample_columns = [c for c in genes.columns if c not in ("gene_id", "gene_name")]
    missing_samples = [c for c in sample_columns if c not in name_map.index]
    if missing_samples:
        raise ValueError(f"TSV samples missing from Excel: {missing_samples[:20]}")

    new_names = [name_map[c] for c in sample_columns]
    dupes = pd.Series(new_names)[pd.Series(new_names).duplicated(keep=False)].unique().tolist()
    if dupes:
        raise ValueError(f"New column names are not unique: {dupes[:10]}")

    # Rename and write: genes stay as rows, no transposing needed
    out = genes.rename(columns=dict(zip(sample_columns, new_names)))
    id_cols = ["gene_id", "gene_name"] if KEEP_GENE_NAME and "gene_name" in out.columns else ["gene_id"]
    out = out[id_cols + new_names]

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_FILE, sep="\t", index=False)

    print(f"Saved: {OUTPUT_FILE.resolve()}")
    print(f"Genes: {len(out)} | Samples: {len(new_names)}")
    print(out.iloc[:5, :len(id_cols) + 4].to_string(index=False))


if __name__ == "__main__":
    main()