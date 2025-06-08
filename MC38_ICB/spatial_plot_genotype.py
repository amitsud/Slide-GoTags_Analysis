# Author: Amit Sud
# Date: 1st May 2025
# Description: Generates spatial plots of Rpl18 mutation status and TCR clonotype identity per sample
#              using annotated single-cell spatial data. Highlights NA, WT, MUT, and specific TCR categories.
# Input:
#   - adata_merged_seurat_filtered: AnnData object with:
#       • obs['x'], obs['y']: spatial coordinates
#       • obs['sample_id']: sample identifiers
#       • obs['Rpl18_mut_tcr']: classification categories: 'Rpl18_WT', 'Rpl18_MUT', 'Rpl18_specific_TCR', etc.
# Output:
#   - Spatial plots per sample (optional: saved as PDF)

import matplotlib.pyplot as plt
import numpy as np

# --- Define color mapping for each category ---
color_mapping = {
    'Rpl18_specific_TCR': '#FF0000',   # Vivid red
    'Rpl18_MUT': '#333399',            # Deep blue
    'Other_clonotype': '#FFCCCC',      # Pale red
    'Rpl18_WT': '#A5CEE3',             # Light blue
    'NA': '#D9D9D9'                    # Light grey for missing values
}

# --- Define plotting order to control layer visibility ---
plot_order = ['NA', 'Rpl18_WT', 'Other_clonotype', 'Rpl18_MUT', 'Rpl18_specific_TCR']

# --- Iterate over each sample to generate individual plots ---
for sample_id in adata_merged_seurat_filtered.obs['sample_id'].unique():
    sample_data = adata_merged_seurat_filtered[
        adata_merged_seurat_filtered.obs['sample_id'] == sample_id
    ]

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal', adjustable='datalim')

    # --- Plot NA values first to ensure background visibility ---
    na_subset = sample_data[sample_data.obs['Rpl18_mut_tcr'].isna()]
    ax.scatter(
        na_subset.obs['x'],
        na_subset.obs['y'],
        c=color_mapping['NA'],
        edgecolor='none',
        s=50,
        alpha=0.9,
        label=f"NA (n={na_subset.shape[0]})"
    )

    # --- Plot remaining categories in defined order ---
    for category in plot_order:
        if category != 'NA' and category in color_mapping:
            subset = sample_data[sample_data.obs['Rpl18_mut_tcr'] == category]
            if subset.shape[0] > 0:
                ax.scatter(
                    subset.obs['x'],
                    subset.obs['y'],
                    c=color_mapping[category],
                    edgecolor='black',
                    s=50,
                    alpha=0.9,
                    label=f"{category} (n={subset.shape[0]})",
                    linewidths=0.5
                )

    # --- Add micron scale bar ---
    micron_to_pixel_ratio = 1
    size_microns = 500
    size_pixels = size_microns * micron_to_pixel_ratio

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    bottom_left_x = x_min + (x_max - x_min) * 0.02
    bottom_left_y = y_min + (y_max - y_min) * 0.02

    ax.plot([bottom_left_x, bottom_left_x + size_pixels], [bottom_left_y, bottom_left_y], 'k-', linewidth=2)
    ax.text(
        bottom_left_x + size_pixels / 2,
        bottom_left_y - (y_max - y_min) * 0.015,
        f'{size_microns} µm',
        ha='center',
        va='top',
        fontsize=10
    )

    ax.plot([bottom_left_x, bottom_left_x], [bottom_left_y, bottom_left_y + size_pixels], 'k-', linewidth=2)
    ax.text(
        bottom_left_x - (x_max - x_min) * 0.01,
        bottom_left_y + size_pixels / 2,
        f'{size_microns} µm',
        ha='right',
        va='center',
        fontsize=10,
        rotation=90
    )

    # --- Customize legend ---
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        fontsize=8,
        title='Rpl18_mut_tcr'
    )

    # --- Plot formatting ---
    ax.set_title(f'Spatial Plot for Sample {sample_id}', fontsize=16)
    ax.set_xlabel('Spatial X', fontsize=12)
    ax.set_ylabel('Spatial Y', fontsize=12)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    # plt.savefig(f"Fig2_Spatial_Plot_Rpl18_mut_tcr_Sample_{sample_id}.pdf", format="pdf", dpi=1200)
    plt.show()
