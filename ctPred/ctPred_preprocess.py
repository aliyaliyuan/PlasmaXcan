#Calculates rank-based percentile from IbGs and adds column, merges with Reference Gene epigenomics csv

import pandas as pd
import numpy as np
import h5py
import os
from pathlib import Path

# ------------------------
# Setup paths
# ------------------------
epi = "/home/aliya/Liver/1111/Gene_epigenomics_v0.csv"
ibg_dir = "/home/aliya/Liver/batch1/IbG/ctrl/"
output_dir = "/home/aliya/Liver/ctPred_output/batch1/"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Get all IbG files
ibg_files = sorted([f for f in os.listdir(ibg_dir) if f.endswith('.tsv')])
print(f"Found {len(ibg_files)} IbG files")


# Track statements
successful_files = 0
skipped_existing = 0
skipped_errors = 0

# ------------------------
# 1. Load metadata 
# ------------------------
meta = pd.read_csv("/home/aliya/Liver/batch1/metadata.csv", dtype=str)
meta = meta.dropna(subset=['TSS_enformer_input'])
meta['chromo'] = meta['TSS_enformer_input'].str.split('_').str[0]
meta['chromo'] = meta['chromo'].astype(str).str.strip()

def fix_chr(val):
    val = val.strip()
    if val.startswith("chr"):
        return val
    else:
        return "chr" + val

meta['chromo'] = meta['chromo'].apply(fix_chr)
meta = meta[['external_gene_name', 'chromo', 'TSS_enformer_input']]
meta = meta.rename(columns={'external_gene_name':'gene_name'})
meta_index = meta.set_index('TSS_enformer_input')

print(f"Loaded metadata with {meta_index.shape[0]} rows")

# ------------------------
# 2. Merge IbG with epi 
# ------------------------

#Convert epi to pandas df
pred_df = pd.read_csv(epi)

current_file = 0


for ibg_file in ibg_files:
    current_file += 1
    ibg_path = os.path.join(ibg_dir, ibg_file)
    ibg_basename = Path(ibg_file).stem  # Remove .tsv extension
        
    # Check if output already exists
    out_filename = f"{ibg_basename}_merged.csv"
    out_path = os.path.join(output_dir, out_filename)
        
    if os.path.exists(out_path):
        print(f"Skipping (already exists): {out_filename}")
        skipped_existing += 1
        continue
        
    print(f"Processing: {ibg_basename}")
        
    try:
        with open(ibg_path) as f:
            header = f.readline().strip().replace('"', '').split()
            
        colnames = ["row_id"] + header
            
        train_expr = pd.read_csv(
            ibg_path,
            sep=r"\s+",
            engine="python",
            names=colnames,
            skiprows=1
        )
            
        train_expr = train_expr.drop(columns=["row_id"])
        train_expr["gene_name"] = train_expr["gene_name"].astype(str).str.strip('"')
            
        ctrls = [c for c in train_expr.columns if c.startswith("ctrl")]
        train_expr[ctrls] = train_expr[ctrls].apply(pd.to_numeric, errors="coerce")
        train_expr["mean_expression"] = train_expr[ctrls].mean(axis=1)
        train_expr = train_expr[["gene_name", "mean_expression"]]
            

        merged_df = pred_df.merge(
            train_expr[['gene_name', 'mean_expression']],
            on='gene_name',
            how='left'
        )
            
        merged_df['rank'] = merged_df['mean_expression'].rank(method='average')
        n = len(merged_df)
        merged_df['percentile'] = ((merged_df['rank'] - 1) / (n - 1) * 100) if n > 1 else 0
        merged_df = merged_df.drop(columns=['mean_expression', 'rank'])
            
        #Re-order columns to expected order of ctPred.py
        feature_cols = sorted(
            [c for c in merged_df.columns if c.startswith("feature_")],
            key=lambda x: int(x.split("_")[1])
        )
            
        merged_df = merged_df[['gene_name', 'chromo'] + feature_cols + ['percentile']]
            
        #saved to csv
        merged_df.to_csv(out_path, index=False)
            
        print(f"  ✓ Saved: {out_filename} (Shape: {merged_df.shape})")
        successful_files += 1
            
    except Exception as e:
        print(f"  ⚠ ERROR processing + {ibg_basename}: {str(e)}")
        skipped_errors += 1
        continue

print(f"\n{'='*60}")
print(f"PROCESSING COMPLETE!")
print(f"{'='*60}")
print(f"Successfully generated: {successful_files} new CSV files")
print(f"Skipped (already exist): {skipped_existing} files")
print(f"Skipped (errors):        {skipped_errors} files")
print(f"\nOutput directory: {output_dir}")
print(f"{'='*60}")
