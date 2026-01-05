#Preprocess with IbG.tsv files (not mean_expression files already)
import pandas as pd
import numpy as np
import h5py

# ------------------------
# 1. Load metadata and extract chromo
# ------------------------
meta = pd.read_csv("/home/aliya/Liver/batch1/metadata.csv", dtype=str)

# Ensure TSS_enformer_input exists
meta = meta.dropna(subset=['TSS_enformer_input'])

# Extract raw chromosome
meta['chromo'] = meta['TSS_enformer_input'].str.split('_').str[0]

# --- FORCE chromo format: chr1, chr2, chrX etc ---
meta['chromo'] = meta['chromo'].astype(str).str.strip()

def fix_chr(val):
    val = val.strip()
    if val.startswith("chr"):
        return val
    else:
        return "chr" + val

meta['chromo'] = meta['chromo'].apply(fix_chr)

# Keep relevant columns
meta = meta[['external_gene_name', 'chromo', 'TSS_enformer_input']]

#Rename 'external_gene_name' to 'gene_name'
meta = meta.rename(columns={'external_gene_name':'gene_name'})

# Index by TSS for matching
meta_index = meta.set_index('TSS_enformer_input')

print(f"Loaded metadata with {meta_index.shape[0]} rows")
print("Example chromo values:", meta['chromo'].unique()[:10])

# ------------------------
# 2. Load H5 predictions and map to gene_name
# ------------------------
h5_path = "/home/aliya/Liver/batch1/h5files/HG00096.h5"

gene_vectors = {}   # gene_name -> 5313-vector
gene_chromo = {}    # gene_name -> chromo

skipped = 0
mapped = 0

with h5py.File(h5_path, "r") as f:
    h5_keys = list(f.keys())
    print(f"Found {len(h5_keys)} prediction entries")

    for key in h5_keys:
        clean_key = key.replace("_predictions", "")

        if clean_key not in meta_index.index:
            skipped += 1
            continue

        gene_entry = meta_index.loc[clean_key]

        if isinstance(gene_entry, pd.DataFrame):
            gene_name = gene_entry['gene_name'].iloc[0]
            chromo = gene_entry['chromo'].iloc[0]
        else:
            gene_name = gene_entry['gene_name']
            chromo = gene_entry['chromo']

        arr = np.array(f[key])

        # Must be (4, 5313)
        if arr.ndim != 2 or arr.shape[1] != 5313:
            skipped += 1
            continue

        vec = arr.mean(axis=0)

        gene_vectors[gene_name] = vec
        gene_chromo[gene_name] = chromo
        mapped += 1

print(f"Mapped predictions for {mapped} genes, skipped {skipped} entries")

# ------------------------
# 3. Convert to DataFrame
# ------------------------

gene_names = list(gene_vectors.keys())
feature_matrix = np.vstack(list(gene_vectors.values()))

# Create feature dataframe explicitly
pred_df = pd.DataFrame(
    feature_matrix,
    columns=[f"feature_{i+1}" for i in range(feature_matrix.shape[1])]
)

# Insert gene_name column
pred_df.insert(0, "gene_name", gene_names)

# Insert chromo column safely
pred_df.insert(1, "chromo", pred_df["gene_name"].map(gene_chromo))

# Force chromo to string again
pred_df["chromo"] = pred_df["chromo"].astype(str)

print("Column layout check:")
print(pred_df.head())
print(pred_df.dtypes[:6])


# ------------------------
# 4. Calculating Mean Expression for Genes
# ------------------------

# ------------------------
# 4. Load IbG expression file (HARD-FIXED parsing)
# ------------------------

# First read just the header to know how many ctrl columns exist
with open("/home/aliya/Liver/batch1/IbG/ctrl/B.tsv") as f:
    header = f.readline().strip().replace('"', '').split()

# header looks like:
# ['gene_name', 'ctrl.01', 'ctrl.02', ...]

# We know the first column in data rows is a numeric row index,
# so we prepend a dummy column name for it
colnames = ["row_id"] + header

train_expr = pd.read_csv(
    "/home/aliya/Liver/batch1/IbG/ctrl/B.tsv",
    sep=r"\s+",
    engine="python",
    names=colnames,
    skiprows=1
)

# Drop row index column
train_expr = train_expr.drop(columns=["row_id"])

# Ensure gene_name is string
train_expr["gene_name"] = train_expr["gene_name"].astype(str).str.strip('"')

# Identify ctrl columns
ctrls = [c for c in train_expr.columns if c.startswith("ctrl")]

# Convert expression columns to numeric
train_expr[ctrls] = train_expr[ctrls].apply(pd.to_numeric, errors="coerce")

# Compute mean expression
train_expr["mean_expression"] = train_expr[ctrls].mean(axis=1)

# Keep only needed columns
train_expr = train_expr[["gene_name", "mean_expression"]]

print("Expression preview:")
print(train_expr.head())
print(train_expr.describe())
print(train_expr.dtypes)

# ------------------------
# 5. Merge predictions with mean expression first
# ------------------------
merged_df = pred_df.merge(
    train_expr[['gene_name', 'mean_expression']],
    on='gene_name',
    how='left'
)

# ------------------------
# 6. Convert mean_expression to rank-based percentiles (only for genes in dataset)
# ------------------------
# Rank the expression values ONLY for genes that are in the merged dataset
merged_df['rank'] = merged_df['mean_expression'].rank(method='average')

# Convert rank to percentile (0-100 scale)
# Percentile = (rank - 1) / (n - 1) * 100
n = len(merged_df)
merged_df['percentile'] = ((merged_df['rank'] - 1) / (n - 1)) if n > 1 else 0

# Drop intermediate columns
merged_df = merged_df.drop(columns=['mean_expression', 'rank'])

print("Percentile conversion summary:")
print(merged_df['percentile'].describe())

# ------------------------
# 7. Reorder columns
# ------------------------
feature_cols = sorted(
    [c for c in merged_df.columns if c.startswith("feature_")],
    key=lambda x: int(x.split("_")[1])
)

merged_df = merged_df[['gene_name', 'chromo'] + feature_cols + ['percentile']]

# ------------------------
# 8. Save output
# ------------------------
out_path = "/home/aliya/Liver/1111/HG00096_B_train.csv"
merged_df.to_csv(out_path, index=False)

print(f"\nSaved final dataset → {out_path}")
print(f"Shape: {merged_df.shape}")
print("Final preview:")
print(merged_df.head())
print("\nFinal column types:")
print(merged_df.dtypes[:5])
print(f"\nPercentile range in final dataset:")
print(f"Min: {merged_df['percentile'].min():.2f}")
print(f"Max: {merged_df['percentile'].max():.2f}")
print(f"Median: {merged_df['percentile'].median():.2f}")
