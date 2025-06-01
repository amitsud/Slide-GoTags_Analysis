# Author: Amit Sud
# Date: 1st May 2025
# Description: Reads a UMAP annotation CSV file and generates a UMAP plot of annotated cell types,
#              mapping detailed cluster labels to broader categories and applying a defined color palette.
#              The output is a high-resolution PDF showing the spatial distribution of key cell types.
# Input:
#   - CSV file with columns: Cluster_Name, umap_1, umap_2
# Output:
#   - PDF file: "Fig1_UMAP_OT1_20k.pdf" saved to output directory

import os
import pandas as pd
import matplotlib.pyplot as plt

# Define output directory and ensure it exists
output_dir = "./"
os.makedirs(output_dir, exist_ok=True)

# Define path to CSV file with annotations and UMAP coordinates
# 🔸 Replace the string below with your actual file path
csv_file_path = "path/to/your/UMAP_annotation_file.csv"

# Load CSV into DataFrame
adata_merged_filtered_umap_seurat = pd.read_csv(csv_file_path)

# Define mapping from Cluster_Name to broader manual cell type annotations
cluster_mapping = {
    "T cells": "T-cells",
    "DC": "Myeloid cells",
    "Fibroblast": "Fibroblasts",
    "Tumor": "Tumor cells",
    "Macrophage": "Myeloid cells",
    "Macrophage proliferating": "Myeloid cells",
    "Myeloid": "Myeloid cells"
}

# Apply mapping to create a new column with simplified cell types
adata_merged_filtered_umap_seurat["manual_celltype_annotation"] = adata_merged_filtered_umap_seurat["Cluster_Name"].map(cluster_mapping)

# Define pastel color palette for each broad cell type
color_palette = {
    "T-cells": '#E32128',         # Vivid Red
    "Tumor cells": '#A5CEE3',     # Light Blue
    "Myeloid cells": '#F57F20',   # Vivid Orange
    "Fibroblasts": '#B1CC70'      # Pastel Green
}

# Define order in which to plot cell types to control layering
plot_order = ["T-cells", "Myeloid cells", "Tumor cells", "Fibroblasts"]

# Calculate total number of cells for summary in legend
total_cells = adata_merged_filtered_umap_seurat.shape[0]

# Initialize UMAP plot
fig, ax = plt.subplots(figsize=(10, 10))

# Plot each cell type according to defined order and colors
for cell_type in plot_order:
    subset = adata_merged_filtered_umap_seurat[
        adata_merged_filtered_umap_seurat["manual_celltype_annotation"] == cell_type
    ]
    if not subset.empty:
        ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            c=color_palette[cell_type],
            edgecolor='black',
            s=50,
            alpha=0.9,
            label=f"{cell_type} (n = {subset.shape[0]})",
            linewidths=0.5
        )

# Add total cell count as first legend entry
handles, labels = ax.get_legend_handles_labels()
handles.insert(0, plt.Line2D([], [], marker='o', color='w', label=f'Total n = {total_cells}', markersize=0))
labels.insert(0, f'Total n = {total_cells}')

# Customize and position legend
legend = ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.025, 1), loc='upper left',
                   fontsize=8, labelspacing=1)
legend.get_frame().set_edgecolor('black')
legend.get_frame().set_linewidth(0)

# Add plot title and clean up axes
ax.set_title("UMAP of Clusters by Cluster_Name_Broad", fontsize=22, pad=20)
ax.set_xlabel("UMAP 1", fontsize=14, labelpad=10)
ax.set_ylabel("UMAP 2", fontsize=14, labelpad=10)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])

# Thicken border of plot box
for spine in ax.spines.values():
    spine.set_linewidth(2)

# Save the plot as a high-resolution PDF
pdf_path = os.path.join(output_dir, "Fig1_UMAP_OT1_20k.pdf")
plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
print(f"Saved plot for sample as {pdf_path}")

# Close the figure to free memory
plt.close(fig)
