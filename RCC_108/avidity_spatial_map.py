# Author: Amit Sud
# Description: This script visualizes spatial maps of T cells grouped by their NSUN5-specific avidity
#              (Strong, Weak, or Other) in spatial transcriptomic samples.
# Input:
#   - AnnData object (adata_merged_filtered_seurat) with:
#       - 'sample_id', 'TRA_cdr3', 'TRB_cdr3', 'x', 'y' in obs
#       - 'NSUN5_avidity_TRA_cdr3' and 'NSUN5_avidity_TRB_cdr3' already assigned
# Output:
#   - PDF spatial plots per sample showing avidity category

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# --- Step 1: Define Strong and Weak Avidity Masks ---
strong_mask = (
    (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'] < 80)) |
    (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'] < 80))
)

weak_mask = (
    ((adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'] >= 80)) &
     (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'] >= 80))) |
    ((adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].isna()) & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'] >= 80)) |
    ((adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].isna()) & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'] >= 80))
)

# --- Step 2: Assign Avidity Category ---
adata_merged_filtered_seurat.obs['Avidity_Category'] = 'Other'
adata_merged_filtered_seurat.obs.loc[strong_mask, 'Avidity_Category'] = 'Strong avidity'
adata_merged_filtered_seurat.obs.loc[weak_mask, 'Avidity_Category'] = 'Weak avidity'

# --- Step 3: Define Colors and Plotting Order ---
avidity_colors = {
    'Strong avidity': '#f76307',   # Red
    'Weak avidity': '#f7cdb2',     # Light Orange
    'Other': '#D3D3D3'             # Light Gray
}
plot_order = ['Other', 'Weak avidity', 'Strong avidity']  # Plot order (back to front)

# --- Step 4: Ensure Spatial Coordinates Exist ---
if 'x' in adata_merged_filtered_seurat.obs and 'y' in adata_merged_filtered_seurat.obs:
    adata_merged_filtered_seurat.obsm['spatial'] = adata_merged_filtered_seurat.obs[['x', 'y']].values
else:
    raise ValueError("Error: 'x' and/or 'y' columns not found in adata_merged_filtered_seurat.obs.")

# --- Step 5: Define Plotting Function ---
def plot_spatial_avidity(adata, sample_id, save_path="."):
    subset = adata[adata.obs['sample_id'] == sample_id]
    plot_data = subset.obs[['x', 'y', 'Avidity_Category']].copy()

    fig, ax = plt.subplots(figsize=(10, 10))
    handles = []

    for category in plot_order:
        color = avidity_colors[category]
        cat_data = plot_data[plot_data['Avidity_Category'] == category]
        scatter = ax.scatter(
            cat_data['x'], cat_data['y'],
            c=color, s=100, alpha=0.9,
            edgecolors='black' if category in ['Strong avidity', 'Weak avidity'] else 'none',
            label=f"{category} (n={cat_data.shape[0]})",
            linewidths=0.5
        )
        handles.append(scatter)

    legend = ax.legend(handles=handles[::-1], bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, frameon=True)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(1.5)

    # Scale bar
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    scale_length_microns = 500
    micron_to_pixel_ratio = 1
    scale_length_pixels = scale_length_microns * micron_to_pixel_ratio
    bottom_left_x = x_min + 0.05 * (x_max - x_min)
    bottom_left_y = y_min + 0.08 * (y_max - y_min)
    ax.plot([bottom_left_x, bottom_left_x + scale_length_pixels], [bottom_left_y, bottom_left_y], 'k-', linewidth=2)
    ax.text(bottom_left_x + scale_length_pixels / 2, bottom_left_y - 0.015 * (y_max - y_min),
            f"{scale_length_microns} µm", ha='center', va='top', fontsize=12)

    # Plot styling
    ax.set_aspect('equal')
    ax.set_title(f'Spatial map of NSUN5 avidity categories\nSample: {sample_id}', fontsize=16, pad=20)
    ax.set_xlabel("Spatial 1", fontsize=14)
    ax.set_ylabel("Spatial 2", fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    # Save as PDF
    pdf_filename = f"{save_path}/spatial_map_avidity_{sample_id}.pdf"
    plt.savefig(pdf_filename, dpi=1200, bbox_inches='tight')
    plt.close(fig)
    print(f"📄 Saved: {pdf_filename}")

# --- Step 6: Generate Plots for All Samples ---
for sample in adata_merged_filtered_seurat.obs['sample_id'].unique():
    plot_spatial_avidity(adata_merged_filtered_seurat, sample)
