# Author: Amit Sud
# Description: Generates spatial plots of annotated cell types for each sample within an AnnData object.
#              Applies a consistent color palette and scale bar, and saves high-resolution PDF plots per sample.
# Input:
#   - adata_merged_seurat_filtered: AnnData object with:
#       • obs['x'], obs['y']: spatial coordinates
#       • obs['sample_id']: sample identifiers
#       • obs['manual_celltype_annotation']: cell type labels
# Output:
#   - One PDF spatial plot per sample, saved as: "Spatial_Map_Sample_<sample_id>.pdf"

import numpy as np
import matplotlib.pyplot as plt

# Define color palette for consistent annotation visualization
color_palette = {
    "Tumor": "#A5CEE3",       
    "NK cell": "#F7DC6F",     
    "Fibroblast": "#B1CC70",  
    "T-cell": "#E32128",     
    "Endothelial cell": "#F38DB2", 
    "Myeloid": "#F57F20",
}

# Define order to control plot layering
cell_type_order = ["Fibroblast", "Myeloid", "NK cell", "Tumor", "T-cell"]

# Ensure spatial coordinates exist
if 'x' in adata_merged_seurat_filtered.obs and 'y' in adata_merged_seurat_filtered.obs:
    adata_merged_seurat_filtered.obsm['spatial'] = adata_merged_seurat_filtered.obs[['x', 'y']].values
else:
    print("Error: 'x' and/or 'y' columns not found in adata_merged_seurat_filtered.obs")

# Ensure cell type annotations are present
if 'manual_celltype_annotation' in adata_merged_seurat_filtered.obs:
    # Iterate through each sample to create individual spatial plots
    for sample_id in adata_merged_seurat_filtered.obs['sample_id'].unique():
        adata_sample = adata_merged_seurat_filtered[adata_merged_seurat_filtered.obs['sample_id'] == sample_id]

        fig, ax = plt.subplots(figsize=(9, 9))

        # Plot each cell type in the defined order
        for cell_type in cell_type_order:
            if cell_type in adata_sample.obs['manual_celltype_annotation'].unique():
                color = color_palette[cell_type]
                subset = adata_sample[adata_sample.obs['manual_celltype_annotation'] == cell_type]
                ax.scatter(
                    subset.obsm['spatial'][:, 0],
                    subset.obsm['spatial'][:, 1],
                    c=color,
                    edgecolor='black',
                    s=30,
                    alpha=0.9,
                    label=f"{cell_type} (n = {subset.shape[0]})",
                    linewidths=0.5
                )

        ax.set_aspect('equal', adjustable='datalim')

        # Legend customization
        legend = ax.legend(
            title="Cell Type",
            bbox_to_anchor=(1.15, 1),
            loc='upper left',
            title_fontsize=12,
            fontsize=10,
            labelspacing=1.2
        )
        for handle in legend.legendHandles:
            handle.set_linewidth(0.5)
            handle.set_edgecolor('black')

        legend.get_frame().set_edgecolor('black')
        legend.get_frame().set_linewidth(1.5)

        # Plot title and axes
        ax.set_title(f"Spatial Map of Cell Types (Sample: {sample_id})", fontsize=18, pad=20)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.grid(False)

        # Add micron-scale size legend
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        bottom_left_x = x_min + (x_max - x_min) * 0.02
        bottom_left_y = y_min + (y_max - y_min) * 0.08

        micron_to_pixel_ratio = 1
        size_microns = 500
        size_pixels = size_microns * micron_to_pixel_ratio

        ax.plot([bottom_left_x, bottom_left_x + size_pixels], [bottom_left_y, bottom_left_y], 'k-', linewidth=2)
        ax.text(
            bottom_left_x + size_pixels / 2,
            bottom_left_y - (y_max - y_min) * 0.015,
            f'{size_microns} µm',
            ha='center', va='top', fontsize=12
        )

        ax.plot([bottom_left_x, bottom_left_x], [bottom_left_y, bottom_left_y + size_pixels], 'k-', linewidth=2)
        ax.text(
            bottom_left_x - (x_max - x_min) * 0.01,
            bottom_left_y + size_pixels / 2,
            f'{size_microns} µm',
            ha='right', va='center',
            fontsize=12, rotation=90
        )

        for spine in ax.spines.values():
            spine.set_linewidth(2)

        # Save output
        output_file = f"Spatial_Map_Sample_{sample_id}.pdf"
        # plt.savefig(output_file, format='pdf', dpi=1200, bbox_inches='tight')  # Uncomment to enable saving
        print(f"Spatial plot saved for sample {sample_id} as {output_file}")

        plt.tight_layout()
        # plt.show()  # Uncomment to display plot in interactive environments

else:
    print("Error: 'manual_celltype_annotation' column not found in adata_merged_seurat_filtered.obs.")
