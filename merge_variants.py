#!/usr/bin/env python3
import pandas as pd
import glob
import os
import re
import argparse
import sys

def natural_key(string_):
    """
    Function for natural sorting (e.g., LGC25-DV14 comes before LGC25-DV137).
    Also ensures samples starting with 'S' are sorted at the end.
    """
    if string_.startswith('S'):
        return [1, string_]
    return [0, [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]]

def main():
    # --- 1. Argument Parser Configuration ---
    parser = argparse.ArgumentParser(description="Merge individual variant annotation files into a unified matrix.")
    
    parser.add_argument("-i", "--input", required=True, 
                        help="Input directory containing *.out files (e.g., ./results/annotated)")
    parser.add_argument("-o", "--output_prefix", default="Yeast_merged", 
                        help="Prefix for the output CSV files (default: Yeast_merged)")
    parser.add_argument("-t", "--threads", type=int, default=1, 
                        help="Number of threads (placeholder for future scaling)")

    args = parser.parse_args()

    # --- 2. Path Validation ---
    input_pattern = os.path.join(args.input, "*.out")
    file_list = sorted(glob.glob(input_pattern))
    
    if not file_list:
        print(f"Error: No files found matching pattern: {input_pattern}")
        sys.exit(1)

    print(f"Found {len(file_list)} files. Starting merge process...")

    # --- 3. Task Definitions ---
    target_columns = {
        "af": "AF",
        "gene": "gene",
        "codon_change": "codon_change",
        "variant_aa": "variant_aa"
    }

    # --- 4. Processing Loops ---
    for col, label in target_columns.items():
        print(f"Processing attribute: {label}...")
        dfs = []
        
        for file in file_list:
            filename = os.path.basename(file)
            # Assuming filename format: SampleName_...
            sample_name = filename.split("_")[0]
            
            try:
                # Read data and ensure position is treated as string
                df = pd.read_csv(file, sep="\t", usecols=["pos", "chrom", col], dtype={'pos': str})
                df = df.rename(columns={col: sample_name})
                dfs.append(df)
            except Exception as e:
                print(f"Warning: Skipping {filename} due to error: {e}")

        if not dfs:
            continue

        # Merge and Pivot
        merged_df = pd.concat(dfs, axis=0)
        merged_df = merged_df.pivot_table(
            index=["pos", "chrom"],
            aggfunc="first"
        ).reset_index()

        # Custom Column Sorting
        samples = [c for c in merged_df.columns if c not in ["pos", "chrom"]]
        sorted_samples = sorted(samples, key=natural_key)
        
        final_columns = ["pos", "chrom"] + sorted_samples
        merged_df = merged_df[final_columns]

        # Fill Missing Values
        if col == "af":
            merged_df = merged_df.fillna(0)
        else:
            merged_df = merged_df.fillna(".")

        # Save Output
        output_name = f"{args.output_prefix}_{label}_matrix.csv"
        merged_df.to_csv(output_name, index=False)
        print(f"Success: Matrix saved to -> {output_name}")

if __name__ == "__main__":
    main()
    print("\n--- All tasks completed successfully ---")