# Author: Amit Sud
# Date: 1st May 2025
# Description: Generates a spatial plot of annotated cell types from a single-cell dataset (`adata_filtered`)
#              by rotating spatial coordinates and coloring each cell type using a custom palette.
#              A custom legend and spatial size scale are added, and the figure is saved as a high-quality PDF.
# Input:
#   - adata_filtered: AnnData object with the following required:
#       • adata_filtered.obs['x'], adata_filtered.obs['y']: spatial coordinates of cells
#       • adata_filtered.obs['manual_celltype_annotation']: assigned cell types per cell
# Output:
#   - A formatted PDF figure named "Fig_1_spatial_plot_formatted.pdf" showing spatial layout by cell type.

# Required packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Define a vivid color palette for specific cell types
color_palette = {
    "T-cells": "#E32128",            # Vivid Red
    "Tumor cells": "#A5CEE3",        # Light Blue
    "Myeloid cells": "#F57F20",      # Vivid Orange
    "Fibroblasts": "#B1CC70",        # Pastel Green
    "Endothelial cells": "#DE4298",  # Vivid Pink
    "Adipocytes": "#7E60A8"          # Deep Purple
}

# Set order for plotting cell types (determines layering and legend order)
cell_type_order = ["Adipocytes", "Endothelial cells", "Fibroblasts", "Myeloid cells", "Tumor cells", "T-cells"]

# Validate spatial coordinate presence
if 'x' in adata_filtered.obs and 'y' in adata_filtered.obs:
    adata_filtered.obsm['spatial'] = adata_filtered.obs[['x', 'y']].values
else:
    print("Error: 'x' and/or 'y' columns not found in adata_filtered.obs")

# Check that cell type annotations are present
if 'manual_celltype_annotation' in adata_filtered.obs:
    # Rotate coordinates by -60 degrees for improved visualization
    angle_rad = np.deg2rad(-60)
    rotation_matrix = np.array([
        [np.cos(angle_rad), -np.sin(angle_rad)],
        [np.sin(angle_rad),  np.cos(angle_rad)]
    ])
    adata_filtered.obsm['spatial_rotated'] = np.dot(adata_filtered.obsm['spatial'], rotation_matrix)

    # Create the figure
    fig, ax = plt.subplots(figsize=(9, 9))

    # Plot each cell type using defined color and order
    for cell_type in cell_type_order:
        if cell_type in adata_filtered.obs['manual_celltype_annotation'].unique():
            color = color_palette[cell_type]
            subset = adata_filtered[adata_filtered.obs['manual_celltype_annotation'] == cell_type]
            ax.scatter(
                subset.obsm['spatial_rotated'][:, 0],
                subset.obsm['spatial_rotated'][:, 1],
                c=color,
                edgecolor='black',
                s=50,
                alpha=0.9,
                linewidths=0.5
            )

    ax.set_aspect('equal', adjustable='datalim')

    # Manually build the legend with cell counts
    handles = []
    for cell_type in cell_type_order:
        if cell_type in adata_filtered.obs['manual_celltype_annotation'].unique():
            color = color_palette[cell_type]
            count = (adata_filtered.obs['manual_celltype_annotation'] == cell_type).sum()
            label = f"{cell_type} (n = {count})"
            handles.append(Patch(facecolor=color, edgecolor='black', label=label, linewidth=0.5))

    legend = ax.legend(
        handles=handles,
        title="Cell Type",
        bbox_to_anchor=(1.15, 1),
        loc='upper left',
        title_fontsize=12,
        fontsize=10,
        labelspacing=1.2
    )
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(1.5)

    # Axis and title formatting
    ax.set_title("Spatial Map of Cell Types", fontsize=18, pad=20)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.grid(False)

    # Calculate bottom-left corner for scale bar
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    bottom_left_x = x_min + 0.02 * (x_max - x_min)
    bottom_left_y = y_min + 0.08 * (y_max - y_min)

    # Define micron scale and convert to pixels
    micron_to_pixel_ratio = 1  # Update this if actual resolution is known
    size_microns = 500
    size_pixels = size_microns * micron_to_pixel_ratio

    # Draw horizontal and vertical scale bars
    ax.plot([bottom_left_x, bottom_left_x + size_pixels],
            [bottom_left_y, bottom_left_y], 'k-', linewidth=2)
    ax.text(bottom_left_x + size_pixels / 2,
            bottom_left_y - 0.015 * (y_max - y_min),
            f'{size_microns} µm', ha='center', va='top', fontsize=12)

    ax.plot([bottom_left_x, bottom_left_x],
            [bottom_left_y, bottom_left_y + size_pixels], 'k-', linewidth=2)
    ax.text(bottom_left_x - 0.01 * (x_max - x_min),
            bottom_left_y + size_pixels / 2,
            f'{size_microns} µm', ha='right', va='center', fontsize=12, rotation=90)

    # Thicken plot borders
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    # Save plot
    output_file = "Fig_1_spatial_plot_formatted.pdf"
    plt.savefig(output_file, format='pdf', bbox_inches='tight')

    plt.tight_layout()
    plt.show()

else:
    print("Error: 'manual_celltype_annotation' column not found in adata_filtered.obs.")
