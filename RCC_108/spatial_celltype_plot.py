# Author: Amit Sud  
# Date: 1st May 2025  
# Description: This script generates spatial scatter plots of annotated cell types 
#              from an AnnData object, grouped by sample. It overlays spatial coordinates, 
#              colors cells by cell type, adds count-based legends, and includes a scale bar.
# Input:
#   - AnnData object with `x`, `y`, `sample_id`, and `manual_celltype_annotation` in `.obs`
# Output:
#   - One spatial PDF plot per sample saved to the specified output directory

import os
import matplotlib.pyplot as plt

# Define output directory for PDFs
output_dir = "path/to/output_directory"
os.makedirs(output_dir, exist_ok=True)

# Define a new color palette for cell types
color_palette = {
    "Tumor": "#A5CEE3",       
    "NK cell": "#F7DC6F",     
    "Fibroblast": "#B1CC70",  
    "T-cell": "#E32128",     
    "Endothelial cell": "#F38DB2", 
    "Myeloid": "#F57F20",
    "B-cell": '#6A0DAD',
}

# Define the plotting order (from back to front)
plot_order = ["Fibroblast", "Endothelial cell", "Fibroblast", "Myeloid", "B-cell", "NK cell", "Tumor", "T-cell"]

# Ensure spatial coordinates are correctly stored in obsm['spatial']
if 'x' in adata_merged_filtered_seurat.obs and 'y' in adata_merged_filtered_seurat.obs:
    coordinates = adata_merged_filtered_seurat.obs[['x', 'y']].values
    adata_merged_filtered_seurat.obsm['spatial'] = coordinates
else:
    print("Error: 'x' and/or 'y' columns not found in adata_merged_filtered_seurat.obs")

# Check if necessary metadata is present
if 'manual_celltype_annotation' in adata_merged_filtered_seurat.obs and 'sample_id' in adata_merged_filtered_seurat.obs:
    unique_samples = adata_merged_filtered_seurat.obs['sample_id'].unique()

    for sample in unique_samples:
        fig, ax = plt.subplots(figsize=(10, 10))
        subset_data = adata_merged_filtered_seurat[adata_merged_filtered_seurat.obs['sample_id'] == sample]
        total_cells = subset_data.shape[0]
        handles = []

        for cell_type in plot_order:
            if cell_type in subset_data.obs['manual_celltype_annotation'].unique():
                subset = subset_data[subset_data.obs['manual_celltype_annotation'] == cell_type]
                color = color_palette[cell_type]
                count = subset.shape[0]
                scatter = ax.scatter(subset.obsm['spatial'][:, 0], subset.obsm['spatial'][:, 1],
                                     c=color, edgecolor='black', s=50, alpha=0.9, 
                                     label=f"{cell_type} (n = {count})", linewidths=0.5)
                handles.append(scatter)
            else:
                handles.append(plt.Line2D([], [], marker='o', color=color_palette[cell_type], linestyle='',
                                          markersize=5, label=f"{cell_type} (n = 0)"))

        # Add total count at top of legend
        handles.insert(0, plt.Line2D([], [], marker='o', color='w', label=f'Total n = {total_cells}', markersize=0))
        legend = ax.legend(handles=handles, bbox_to_anchor=(1.025, 1), loc='upper left', fontsize=8, labelspacing=1)
        legend.get_frame().set_edgecolor('black')
        legend.get_frame().set_linewidth(0)

        ax.set_title(f"Spatial Map of Cell Types (Sample: {sample})", fontsize=18, pad=20)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.grid(False)
        ax.set_aspect('equal', adjustable='datalim')

        # Scale bar
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        bottom_left_x = x_min + (x_max - x_min) * 0.02
        bottom_left_y = y_min + (y_max - y_min) * 0.08
        micron_to_pixel_ratio = 1
        size_microns = 500
        size_pixels = size_microns * micron_to_pixel_ratio

        ax.plot([bottom_left_x, bottom_left_x + size_pixels], [bottom_left_y, bottom_left_y], 'k-', linewidth=2)
        ax.text(bottom_left_x + size_pixels / 2, bottom_left_y - (y_max - y_min) * 0.015, f'{size_microns} µm',
                ha='center', va='top', fontsize=12)

        ax.plot([bottom_left_x, bottom_left_x], [bottom_left_y, bottom_left_y + size_pixels], 'k-', linewidth=2)
        ax.text(bottom_left_x - (x_max - x_min) * 0.01, bottom_left_y + size_pixels / 2, f'{size_microns} µm',
                ha='right', va='center', fontsize=12, rotation=90)

        for spine in ax.spines.values():
            spine.set_linewidth(2)

        pdf_path = os.path.join(output_dir, f"spatial_plot_sample_{sample}.pdf")
        # Uncomment below to save the plots
        # plt.savefig(pdf_path, format='pdf', dpi=1200, bbox_inches='tight')
        plt.show()
        plt.close(fig)
        print(f"Saved plot for sample {sample} as {pdf_path}")

else:
    print("Required columns 'manual_celltype_annotation' or 'sample_id' not found in adata_merged_filtered_seurat.obs.")

