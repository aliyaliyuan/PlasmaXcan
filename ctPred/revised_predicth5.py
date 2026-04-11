#!/usr/bin/env python3

#Matches by genomic overlap

import torch
import pandas as pd
import h5py
import numpy as np
import os
import sys
import time
from pathlib import Path
import argparse

def parse_region(region_str):
    """Parse region string to chr, start, end"""
    # Remove _predictions suffix
    region_str = region_str.replace('_predictions', '')
    
    # Format: chr10_100009946_100009948
    parts = region_str.split('_')
    if len(parts) == 3:
        chr_name = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        return chr_name, start, end
    return None, None, None

def regions_overlap(chr1, start1, end1, chr2, start2, end2):
    """Check if two regions overlap"""
    if chr1 != chr2:
        return False
    # Check overlap
    return not (end1 < start2 or end2 < start1)

def find_gene_for_region(region_str, ref_df, tss_col='TSS_enformer_input', chr_col='chromosome', tss_pos_col='TSS'):
    """Find matching gene by checking if region contains TSS"""
    chr_h5, start_h5, end_h5 = parse_region(region_str)
    
    if chr_h5 is None:
        return None
    
    # Filter to same chromosome
    chr_matches = ref_df[ref_df[chr_col] == chr_h5]
    
    # Find genes where TSS falls within H5 region
    for idx, row in chr_matches.iterrows():
        tss = int(row[tss_pos_col])
        
        # Check if TSS is within the H5 region
        if start_h5 <= tss <= end_h5:
            return row['ensembl_gene_id']
    
    return None

def print_progress(current, total, start_time, prefix=""):
    """Simple progress display"""
    elapsed = time.time() - start_time
    percent = 100 * (current / total) if total > 0 else 0
    bar_length = 40
    filled_length = int(bar_length * current // total)
    bar = '█' * filled_length + '░' * (bar_length - filled_length)
    
    sys.stdout.write(f'\r{prefix} |{bar}| {current}/{total} ({percent:.1f}%) {elapsed:.1f}s')
    sys.stdout.flush()
    
    if current == total:
        print()

def main():
    from ctPred_utils import ctPred
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Predict gene expression from H5 files')
    parser.add_argument('--model_path', required=True, help='Path to trained model (.pt file)')
    parser.add_argument('--input_folder', required=True, help='Folder containing .h5 files')
    parser.add_argument('--output_path', required=True, help='Output CSV file path')
    parser.add_argument('--ref_csv', required=True, help='Reference CSV file with gene mappings')
    parser.add_argument('--device', default='cpu', help='Device to use (cpu or cuda)')
    args = parser.parse_args()
    
    # Load model
    device = torch.device(args.device)
    model = ctPred().to(device)
    model.load_state_dict(torch.load(args.model_path, map_location=device))
    model.eval()
    
    # Load reference CSV
    print("Loading reference CSV...")
    ref_df = pd.read_csv(args.ref_csv)
    
    # Clean up ensembl_gene_id column (remove version numbers if present)
    ref_df['ensembl_gene_id_clean'] = ref_df['ensembl_gene_id'].str.split('.').str[0]
    
    print(f"Reference genes: {len(ref_df)}")
    print(f"Chromosomes: {ref_df['chromosome'].unique()[:5]}...")
    
    # Get all H5 files
    h5_files = list(Path(args.input_folder).glob("*.h5"))
    
    if not h5_files:
        print(f"Error: No .h5 files found in {args.input_folder}")
        sys.exit(1)
    
    # Process each file
    results = {}
    total_genes = 0
    start_time = time.time()
    
    print(f"\nProcessing {len(h5_files)} samples...")
    
    for i, h5_file in enumerate(h5_files, 1):
        print_progress(i, len(h5_files), start_time, "Samples")
        
        with h5py.File(h5_file, 'r') as f:
            preds = {}
            gene_keys = list(f.keys())
            
            for gene_key in gene_keys:
                data = f[gene_key][:]
                features = np.mean(data, axis=0)
                
                X = torch.tensor(features, dtype=torch.float32).to(device)
                with torch.no_grad():
                    pred = model(X.unsqueeze(0)).item()
                
                # Store with original key (will map later)
                preds[gene_key] = pred
            
            results[h5_file.stem] = preds
            total_genes += len(preds)
    
    print(f"\nProcessed {total_genes} predictions")
    
    # Create DataFrame
    df = pd.DataFrame.from_dict(results, orient='index').T
    
    # Map regions to gene IDs
    print("\nMapping regions to genes...")
    region_to_gene = {}
    matched = 0
    unmatched = 0
    
    total_regions = len(df.index)
    mapping_start = time.time()
    
    for i, region in enumerate(df.index, 1):
        if i % 100 == 0 or i == total_regions:
            print_progress(i, total_regions, mapping_start, "Mapping")
        
        gene_id = find_gene_for_region(region, ref_df)
        
        if gene_id:
            # Remove version number from gene ID if present
            gene_id_clean = gene_id.split('.')[0]
            region_to_gene[region] = gene_id_clean
            matched += 1
        else:
            unmatched += 1
    
    print(f"\n  Matched: {matched}")
    print(f"  Unmatched: {unmatched}")
    
    # Filter DataFrame to only matched regions
    if region_to_gene:
        valid_regions = list(region_to_gene.keys())
        df_matched = df.loc[valid_regions].copy()
        
        # Map to gene IDs
        df_matched.index = df_matched.index.map(region_to_gene)
        
        # Handle duplicates (multiple regions mapping to same gene - keep first)
        df_matched = df_matched[~df_matched.index.duplicated(keep='first')]
        
        df_matched.index.name = 'NAME'
        df_matched = df_matched.reset_index()
        
        # Save
        df_matched.to_csv(args.output_path, index=False)
        
        print(f"\n✓ Completed in {time.time() - start_time:.1f}s")
        print(f"  Output: {args.output_path}")
        print(f"  Unique genes: {len(df_matched)}")
        print(f"  Samples: {df.shape[1]}")
    else:
        print(f"\n✗ Error: No matches found")
        sys.exit(1)

if __name__ == "__main__":
    main()
