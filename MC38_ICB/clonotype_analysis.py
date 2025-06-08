# Author: Amit Sud
# Date: 1st May 2025
# Description: Generates bar plots of TRA and TRB clonotype frequencies from T cells across all samples,
#              and individually per sample. Highlights specific Rpl18-associated clonotypes using a defined color scheme.
# Input:
#   - adata_merged_seurat_filtered: AnnData object with:
#       • obs['TRA_cdr3'], obs['TRB_cdr3']: clonotype CDR3 sequences
#       • obs['manual_celltype_annotation']: must include 'T-cell'
#       • obs['sample_id']: sample identifiers
#       • Mustafa_TRA, Mustafa_TRB, Bulk_TRA, Bulk_TRB: sets of CDR3 sequences for highlighting
# Output:
#   - Bar plots saved as PDF (optional) and displayed using `matplotlib`

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

# --- Parameters ---
bar_width = 0.8
padding_factor = 0.6
dpi_value = 1200
output_dir = "plots"
os.makedirs(output_dir, exist_ok=True)

# --- Subset for T cells only ---
adata_t_cells = adata_merged_seurat_filtered[
    adata_merged_seurat_filtered.obs['manual_celltype_annotation'] == 'T-cell'
]

# === 1. PLOT FOR ALL SAMPLES COMBINED ===
TRA_counts_all = adata_t_cells.obs['TRA_cdr3'].value_counts()
TRB_counts_all = adata_t_cells.obs['TRB_cdr3'].value_counts()

cdr3_tra_all = pd.DataFrame({'TRA': TRA_counts_all})
cdr3_trb_all = pd.DataFrame({'TRB': TRB_counts_all})

if not cdr3_tra_all.empty and not cdr3_trb_all.empty:
    # Color by clonotype group membership
    bar_colors_tra = ['#FF0000' if c in Bulk_TRA else '#FFCCCC' for c in cdr3_tra_all.index]
    bar_colors_trb = ['#FF0000' if c in Bulk_TRB else '#FFCCCC' for c in cdr3_trb_all.index]

    fig_width_tra = len(cdr3_tra_all) * bar_width * padding_factor
    fig_width_trb = len(cdr3_trb_all) * bar_width * padding_factor
    fig = plt.figure(figsize=(fig_width_tra + fig_width_trb + 2, 8))
    gs = fig.add_gridspec(1, 2, width_ratios=[fig_width_tra, fig_width_trb], wspace=0.3)

    # --- Plot TRA ---
    ax1 = fig.add_subplot(gs[0])
    ax1.bar(range(len(cdr3_tra_all)), cdr3_tra_all['TRA'], color=bar_colors_tra, edgecolor='black', width=bar_width)
    ax1.set_title('All Samples - TRA Clonotypes', fontsize=18)
    ax1.set_xlabel('TRA Clonotype', fontsize=14)
    ax1.set_ylabel('Number of T cells', fontsize=14)
    ax1.set_xticks(range(len(cdr3_tra_all)))
    ax1.set_xticklabels(cdr3_tra_all.index, rotation=90, fontsize=10)
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax1.set_xlim(-0.5, len(cdr3_tra_all) - 0.5)
    ax1.grid(False)

    # --- Plot TRB ---
    ax2 = fig.add_subplot(gs[1])
    ax2.bar(range(len(cdr3_trb_all)), cdr3_trb_all['TRB'], color=bar_colors_trb, edgecolor='black', width=bar_width)
    ax2.set_title('All Samples - TRB Clonotypes', fontsize=18)
    ax2.set_xlabel('TRB Clonotype', fontsize=14)
    ax2.set_ylabel('Number of T cells', fontsize=14)
    ax2.set_xticks(range(len(cdr3_trb_all)))
    ax2.set_xticklabels(cdr3_trb_all.index, rotation=90, fontsize=10)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax2.set_xlim(-0.5, len(cdr3_trb_all) - 0.5)
    ax2.grid(False)

    # Equalize y-axis
    y_max = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_ylim(0, y_max)
    ax2.set_ylim(0, y_max)

    # Legend
    legend_patches = [
        Patch(facecolor='#FF0000', edgecolor='black', label='Rpl18 specific clonotype', linewidth=1),
        Patch(facecolor='#FFCCCC', edgecolor='black', label='Other Clonotypes', linewidth=1)
    ]
    fig.legend(handles=legend_patches, title='Clonotype Highlights', fontsize=12, title_fontsize=12,
               loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.05), frameon=True)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    # plt.savefig(f"{output_dir}/TRA_TRB_Clonotypes_All_Samples.pdf", format="pdf", dpi=dpi_value, bbox_inches="tight")
    print(f"Saved TRA_TRB_Clonotypes_All_Samples.pdf")
    plt.show()

# === 2. PLOT PER INDIVIDUAL SAMPLE ===
for sample_id in adata_t_cells.obs['sample_id'].unique():
    adata_sample = adata_t_cells[adata_t_cells.obs['sample_id'] == sample_id]

    TRA_counts = adata_sample.obs['TRA_cdr3'].value_counts()
    TRB_counts = adata_sample.obs['TRB_cdr3'].value_counts()

    cdr3_tra = pd.DataFrame({'TRA': TRA_counts})
    cdr3_trb = pd.DataFrame({'TRB': TRB_counts})

    if cdr3_tra.empty or cdr3_trb.empty:
        print(f"Skipping sample {sample_id}: No TRA or TRB clonotypes.")
        continue

    bar_colors_tra = ['#FF0000' if c in Bulk_TRA else '#FFCCCC' for c in cdr3_tra.index]
    bar_colors_trb = ['#FF0000' if c in Bulk_TRB else '#FFCCCC' for c in cdr3_trb.index]

    fig_width_tra = len(cdr3_tra) * bar_width * padding_factor
    fig_width_trb = len(cdr3_trb) * bar_width * padding_factor
    fig = plt.figure(figsize=(fig_width_tra + fig_width_trb + 2, 8))
    gs = fig.add_gridspec(1, 2, width_ratios=[fig_width_tra, fig_width_trb], wspace=0.3)

    # --- Plot TRA ---
    ax1 = fig.add_subplot(gs[0])
    ax1.bar(range(len(cdr3_tra)), cdr3_tra['TRA'], color=bar_colors_tra, edgecolor='black', width=bar_width)
    ax1.set_title(f'{sample_id} - TRA Clonotypes', fontsize=18)
    ax1.set_xlabel('TRA Clonotype', fontsize=14)
    ax1.set_ylabel('Number of T cells', fontsize=14)
    ax1.set_xticks(range(len(cdr3_tra)))
    ax1.set_xticklabels(cdr3_tra.index, rotation=90, fontsize=10)
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax1.set_xlim(-0.5, len(cdr3_tra) - 0.5)
    ax1.grid(False)

    # --- Plot TRB ---
    ax2 = fig.add_subplot(gs[1])
    ax2.bar(range(len(cdr3_trb)), cdr3_trb['TRB'], color=bar_colors_trb, edgecolor='black', width=bar_width)
    ax2.set_title(f'{sample_id} - TRB Clonotypes', fontsize=18)
    ax2.set_xlabel('TRB Clonotype', fontsize=14)
    ax2.set_ylabel('Number of T cells', fontsize=14)
    ax2.set_xticks(range(len(cdr3_trb)))
    ax2.set_xticklabels(cdr3_trb.index, rotation=90, fontsize=10)
    ax2.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax2.set_xlim(-0.5, len(cdr3_trb) - 0.5)
    ax2.grid(False)

    y_max = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_ylim(0, y_max)
    ax2.set_ylim(0, y_max)

    fig.legend(handles=legend_patches, title='Clonotype Highlights', fontsize=12, title_fontsize=12,
               loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.05), frameon=True)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    # plt.savefig(f"{output_dir}/TRA_TRB_Clonotypes_{sample_id}.pdf", format="pdf", dpi=dpi_value, bbox_inches="tight")
    print(f"Saved TRA_TRB_Clonotypes_{sample_id}.pdf")
    # plt.show()
