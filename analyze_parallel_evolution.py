#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re
import argparse
import sys

# =================================================================
# Helper Functions for Data Cleaning and Classification
# =================================================================

def extract_first_valid(row):
    """Extract the first non-null/non-dot value from a matrix row."""
    for v in row[2:]:
        if pd.notna(v) and v != ".":
            return v
    return "NA"

def clean_codon(c):
    """Clean codon change strings."""
    if c == "NA": return c
    c = str(c).replace("c.", "")
    c = re.sub(r"^-", "", c)
    return c

def clean_aa(a):
    """Translate 3-letter AA codes to 1-letter codes."""
    if a == "NA": return a
    a = str(a).replace("p.", "")
    
    aa_map = {
        "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
        "Glu":"E","Gln":"Q","Gly":"G","His":"H","Ile":"I",
        "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
        "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
        "Ter":"*","Stop":"*"
    }
    
    # Match pattern: (Ref)(Pos)(Alt)
    m = re.match(r"([A-Za-z]+)(\d+)([A-Za-z\*]+)", a)
    if m:
        ref, pos, alt = m.group(1), m.group(2), m.group(3)
        ref1 = aa_map.get(ref, ref)
        alt1 = aa_map.get(alt, alt)
        return f"{ref1}{pos}{alt1}"
    return a

def classify_mutation(aa):
    """Classify mutation as synonymous, nonsynonymous, or nonsense."""
    if aa == "NA": return "unknown"
    types = []
    for a in str(aa).split(","):
        try:
            ref, alt = a[0], a[-1]
            if alt == "*": types.append("nonsense")
            elif ref == alt: types.append("synonymous")
            else: types.append("nonsynonymous")
        except:
            types.append("unknown")
    return ",".join(sorted(set(types)))

# =================================================================
# Main Analysis Function
# =================================================================

def main():
    parser = argparse.ArgumentParser(description="Analyze parallel evolution and generate rankings/plots.")
    
    # Input arguments
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing candidate mutation TSVs")
    parser.add_argument("-m", "--matrix_dir", required=True, help="Directory containing Codon/AA matrices")
    parser.add_argument("-o", "--output_dir", default="parallel_evolution_results", help="Output directory")
    parser.add_argument("--top_n", type=int, default=20, help="Number of top genes to plot (default: 20)")
    parser.add_argument("--af_col", default="T45", help="Column name for Allele Frequency (default: T45)")

    args = parser.parse_args()

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Load Annotation Matrices
    print("Loading annotation matrices...")
    codon_files = glob.glob(os.path.join(args.matrix_dir, "*_codon_change_matrix.csv"))
    aa_files = glob.glob(os.path.join(args.matrix_dir, "*_variant_aa_matrix.csv"))

    if not codon_files or not aa_files:
        print("Error: Missing codon or AA matrices in matrix directory.")
        sys.exit(1)

    codon_map = {}
    for f in codon_files:
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            codon_map[(row["chrom"], row["pos"])] = clean_codon(extract_first_valid(row))

    aa_map = {}
    for f in aa_files:
        df = pd.read_csv(f)
        for _, row in df.iterrows():
            aa_map[(row["chrom"], row["pos"])] = clean_aa(extract_first_valid(row))

    # 2. Load Candidate Mutations
    candidate_files = glob.glob(os.path.join(args.input_dir, "*_candidate_mutations.tsv"))
    if not candidate_files:
        print(f"Error: No candidate mutation files found in {args.input_dir}")
        sys.exit(1)

    dfs = []
    for f in candidate_files:
        df = pd.read_csv(f, sep="\t")
        group = os.path.basename(f).split("_")[0]
        df["Founder"] = group
        dfs.append(df)
    
    data = pd.concat(dfs, ignore_index=True)

    # 3. Add Annotation Columns
    print("Integrating annotations...")
    data["Codon_change"] = [codon_map.get((r["CHR"], r["POS"]), "NA") for _, r in data.iterrows()]
    data["AA_change"] = [aa_map.get((r["CHR"], r["POS"]), "NA") for _, r in data.iterrows()]
    data["Mutation_type"] = data["AA_change"].apply(classify_mutation)

    # 4. Save Detailed SNP Table
    data = data.sort_values(["Gene", "Founder", "Lineage"])
    detailed_file = os.path.join(args.output_dir, "detailed_parallel_evolution_SNPs.tsv")
    data.to_csv(detailed_file, sep="\t", index=False)

    # 5. Generate Gene-level Summary
    print("Generating gene-level summary...")
    summary = data.groupby("Gene").agg(
        Total_SNPs=("Gene", "count"),
        Replicate_hits=("Lineage", "nunique"),
        Founder_hits=("Founder", "nunique"),
        Mean_AF=(args.af_col, "mean")
    ).reset_index()

    # Selection score calculation
    summary["Selection_score"] = (summary["Founder_hits"] * 3 + summary["Replicate_hits"] * 2 + summary["Total_SNPs"])
    summary = summary.sort_values("Selection_score", ascending=False)
    summary.to_csv(os.path.join(args.output_dir, "gene_selection_ranking.tsv"), sep="\t", index=False)

    # 6. Generate Amino Acid-level Summary
    print("Generating AA-level summary...")
    aa_snps = data[data["AA_change"] != "NA"].copy()
    aa_stats = aa_snps.groupby(["Gene", "AA_change"]).agg(
        SNP_count=("AA_change", "count"),
        Founder_hits=("Founder", "nunique"),
        Replicate_hits=("Lineage", "nunique"),
        Mean_AF=(args.af_col, "mean")
    ).reset_index()
    aa_stats["Selection_score"] = (aa_stats["Founder_hits"] * 3 + aa_stats["Replicate_hits"] * 2 + aa_stats["SNP_count"])
    aa_stats = aa_stats.sort_values("Selection_score", ascending=False)
    aa_stats.to_csv(os.path.join(args.output_dir, "aa_parallel_evolution_ranking.tsv"), sep="\t", index=False)

    # 7. Visualization
    print(f"Generating plots for top {args.top_n} genes...")
    top_genes = summary.head(args.top_n)

    # Plot Founder Hits
    plt.figure(figsize=(10, 8))
    plt.barh(top_genes["Gene"][::-1], top_genes["Founder_hits"][::-1], color='skyblue')
    plt.xlabel("Number of Founder Groups")
    plt.title("Parallel Evolution Across Founders")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "plot_founder_parallelism.pdf"))

    # Plot Replicate Hits
    plt.figure(figsize=(10, 8))
    plt.barh(top_genes["Gene"][::-1], top_genes["Replicate_hits"][::-1], color='salmon')
    plt.xlabel("Number of Replicate Lineages")
    plt.title("Parallel Evolution Across Replicates")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "plot_replicate_parallelism.pdf"))

    print(f"\n--- Analysis Complete ---")
    print(f"Results saved in: {args.output_dir}")

if __name__ == "__main__":
    main()