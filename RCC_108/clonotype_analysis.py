# Author: Amit Sud
# Description: This script visualizes TRA and TRB clonotype frequencies from an AnnData object.
#              Clonotypes are grouped and colored by antigen specificity (e.g., NSUN5, Tumor).
#              For each sample, separate bar plots are generated for TRA and TRB CDR3 counts,
#              with consistent axis scaling, color-coded highlights, and saved in high resolution.
# Input:
#   - AnnData object with `TRA_cdr3`, `TRB_cdr3`, `manual_celltype_annotation`, and `sample_id` in `.obs`
#   - A DataFrame (`rcc_tcr_specificity`) specifying clonotype antigen specificity
# Output:
#   - Per-sample PDF plots saved to a specified output directory

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path

# Parameters
bar_width = 0.8  # Fixed bar width
padding_factor = 0.6  # Adjusts space allocated for each bar
dpi_value = 1200  # High-resolution output
output_dir = Path("path/to/output_plots")  # Generic path to save plots
output_dir.mkdir(exist_ok=True)

# Helper function to get color based on prioritized specificity
def get_color(cdr3_seq, column, color_priority):
    for specificity_df, color in color_priority:
        if cdr3_seq in specificity_df[column].values:
            return color
    return '#FFCCCC'  # Default color for other clonotypes

# Define subsets for each specificity in priority order
nsun5_specificity = rcc_tcr_specificity[rcc_tcr_specificity['Specificity'].str.contains("NSUN5", na=False)]
tumor_specificity = rcc_tcr_specificity[rcc_tcr_specificity['Specificity'] == "Tumor"]
color_priority = [
    (nsun5_specificity, '#FF0000'),    # NSUN5 = red
    (tumor_specificity, '#F8C47D')     # Tumor = light orange
]

# Legend entries
legend_patches = [
    Patch(facecolor='#FF0000', edgecolor='black', label='NSUN5 specific clonotype'),
    Patch(facecolor='#F8C47D', edgecolor='black', label='Tumor specific'),
    Patch(facecolor='#FFCCCC', edgecolor='black', label='Other clonotypes'),
]

# Filter AnnData for T-cells
adata_t_cells = adata_merged_filtered_seurat[
    adata_merged_filtered_seurat.obs['manual_celltype_annotation'] == 'T-cell'
]

# Combined counts
combined_tra_counts = None
combined_trb_counts = None

# Function to calculate standardized figure width
def calculate_figure_width(num_bars, bar_width, padding_factor):
    return num_bars * bar_width * padding_factor

# Process each sample
for sample_id in adata_t_cells.obs['sample_id'].unique():
    adata_sample = adata_t_cells[adata_t_cells.obs['sample_id'] == sample_id]
    tra_counts = adata_sample.obs['TRA_cdr3'].value_counts()
    trb_counts = adata_sample.obs['TRB_cdr3'].value_counts()

    # Print counts
    print(f"Sample {sample_id}: Unique TRA_cdr3 = {adata_sample.obs['TRA_cdr3'].nunique()}, "
          f"Unique TRB_cdr3 = {adata_sample.obs['TRB_cdr3'].nunique()}")

    # Aggregate
    combined_tra_counts = tra_counts if combined_tra_counts is None else combined_tra_counts.add(tra_counts, fill_value=0)
    combined_trb_counts = trb_counts if combined_trb_counts is None else combined_trb_counts.add(trb_counts, fill_value=0)

    # Color assignment
    tra_colors = [get_color(cdr3, 'CDR3A', color_priority) for cdr3 in tra_counts.index]
    trb_colors = [get_color(cdr3, 'CDR3B', color_priority) for cdr3 in trb_counts.index]

    max_num_bars = max(len(tra_counts), len(trb_counts))
    fig_width = calculate_figure_width(max_num_bars, bar_width, padding_factor)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, 8), gridspec_kw={'wspace': 0.4})

    # TRA plot
    ax1.bar(range(len(tra_counts)), tra_counts.values, color=tra_colors, edgecolor='black', width=bar_width)
    ax1.set_title(f'{sample_id} - TRA Clonotypes', fontsize=16)
    ax1.set_xlabel('TRA Clonotypes', fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_xticks(range(len(tra_counts)))
    ax1.set_xticklabels(tra_counts.index, rotation=90, fontsize=10)
    ax1.set_xlim(-0.5, max_num_bars - 0.5)
    ax1.grid(False)

    # TRB plot
    ax2.bar(range(len(trb_counts)), trb_counts.values, color=trb_colors, edgecolor='black', width=bar_width)
    ax2.set_title(f'{sample_id} - TRB Clonotypes', fontsize=16)
    ax2.set_xlabel('TRB Clonotypes', fontsize=12)
    ax2.set_ylabel('Count', fontsize=12)
    ax2.set_xticks(range(len(trb_counts)))
    ax2.set_xticklabels(trb_counts.index, rotation=90, fontsize=10)
    ax2.set_xlim(-0.5, max_num_bars - 0.5)
    ax2.grid(False)

    # Consistent y-axis scaling
    max_y = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_ylim(0, max_y)
    ax2.set_ylim(0, max_y)

    # Legend
    fig.legend(handles=legend_patches, title='Clonotype Highlights', loc='upper center', ncol=3, frameon=True)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save plot
    plt.savefig(output_dir / f'TRA_TRB_Clonotypes_{sample_id}.pdf', dpi=dpi_value, bbox_inches='tight')
    plt.show()
    plt.close()

# Global counts
unique_tra_all = adata_t_cells.obs['TRA_cdr3'].dropna().nunique()
unique_trb_all = adata_t_cells.obs['TRB_cdr3'].dropna().nunique()
print(f"Across All Samples: Unique TRA_cdr3 = {unique_tra_all}, Unique TRB_cdr3 = {unique_trb_all}")
print("Plots saved in:", output_dir)
