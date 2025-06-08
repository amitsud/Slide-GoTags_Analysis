# Author: Amit Sud
# Date: 1st May 2025
# Description: This script summarizes cell type distributions across NSUN5 mutant and wild-type genotypes.
#              It generates separate histograms for each genotype, highlighting tumor cell counts.
# Input:
#   - AnnData object with `NSUN5_Genotype_Threshold_Final` and `manual_celltype_annotation` in `.obs`
# Output:
#   - Two high-resolution PDF histograms (one each for NSUN5_MUT and NSUN5_WT)

import matplotlib.pyplot as plt
import pandas as pd

# Group and count cells by genotype and cell type
cell_count_df = adata_merged_filtered_seurat.obs.groupby(
    ['NSUN5_Genotype_Threshold_Final', 'manual_celltype_annotation']
).size().reset_index(name='Cell Count')

# Pivot to wide format: genotypes as rows, cell types as columns
cell_count_df = cell_count_df.pivot(
    index='NSUN5_Genotype_Threshold_Final',
    columns='manual_celltype_annotation',
    values='Cell Count'
).fillna(0).astype(int)

# Define colors
mut_color = '#333399'  # NSUN5_MUT tumor = dark pastel blue
wt_color = '#A5CEE3'   # NSUN5_WT tumor = light pastel blue

# Output directory (generic placeholder)
output_dir = Path("path/to/output_plots")
output_dir.mkdir(exist_ok=True)

# Function to plot histogram for a genotype row
def plot_genotype_histogram(genotype_label, color, output_filename):
    if genotype_label in cell_count_df.index:
        data = cell_count_df.loc[genotype_label]
        fig, ax = plt.subplots(figsize=(10, 6))

        for cluster, count in data.items():
            bar_color = color if cluster == 'Tumor' else 'white'
            ax.bar(cluster, count, color=bar_color, edgecolor='black')

        ax.set_title(f"Histogram of Cell Counts for {genotype_label}", fontsize=14)
        ax.set_xlabel("Cell Type", fontsize=12)
        ax.set_ylabel("Cell Count", fontsize=12)
        ax.tick_params(axis='x', rotation=90)
        ax.grid(False)

        # Save plot
        plt.savefig(output_dir / output_filename, dpi=1200, bbox_inches='tight')
        plt.show()
        plt.close(fig)
        print(f"Histogram saved: {output_filename}")

# Plot and save for both genotypes
plot_genotype_histogram('NSUN5_MUT', mut_color, "NSUN5_Cell_Count_Histogram_MUT.pdf")
plot_genotype_histogram('NSUN5_WT', wt_color, "NSUN5_Cell_Count_Histogram_WT.pdf")
