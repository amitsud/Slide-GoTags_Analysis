# Author: Amit Sud
# Date: 1st May 2025
# Description: This script generates spatial gene expression plots for specified genes and cell populations,
#              overlaid with a KDE-defined region based on selected categories of cells. It creates contour
#              regions using Gaussian KDE and visualizes gene expression with colorbars and scale bar.
# Input:
#   - AnnData object (adata) with:
#       - obs columns: 'sample_id', 'x', 'y', 'NSUN5_MUT_TCR', 'Merged_Cluster_Name_scevan'
#       - expression matrix containing relevant gene names
#   - KDE categories for defining a region
#   - Gene-population pairs for visualization
# Output:
#   - PDF plot showing spatial gene expression with KDE overlay and color legends

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
from shapely.geometry import MultiPoint
import matplotlib.path as mpath
from mpl_toolkits.axes_grid1 import make_axes_locatable

def generate_kde_contour(adata, sample_id, categories, grid_size=250, percentile=90):
    kde_data = adata[
        (adata.obs['sample_id'] == sample_id) &
        (adata.obs['NSUN5_MUT_TCR'].isin(categories))
    ]
    kde_coordinates = kde_data.obs[['x', 'y']].values

    if kde_coordinates.shape[0] == 0:
        print("❌ No data for KDE generation.")
        return None, None, None

    kde = gaussian_kde(kde_coordinates.T, bw_method='scott')
    x_min, x_max = kde_coordinates[:, 0].min(), kde_coordinates[:, 0].max()
    y_min, y_max = kde_coordinates[:, 1].min(), kde_coordinates[:, 1].max()

    x_grid, y_grid = np.meshgrid(
        np.linspace(x_min, x_max, grid_size),
        np.linspace(y_min, y_max, grid_size)
    )
    density = kde(np.vstack([x_grid.ravel(), y_grid.ravel()])).reshape(x_grid.shape)
    threshold = np.percentile(density, percentile)

    contour_set = plt.contour(x_grid, y_grid, density, levels=[threshold])
    plt.close()

    largest_contour = max(contour_set.allsegs[0], key=len, default=None)
    if largest_contour is None:
        print("❌ No contour found.")
        return None, None, None

    centroid = MultiPoint(largest_contour).centroid
    path = mpath.Path(largest_contour)
    return path, centroid.x, centroid.y


def plot_gene_expression_with_kde(adata, sample_id, gene_population_pairs, kde_path=None, save_path="."):
    subset = adata[adata.obs['sample_id'] == sample_id]
    fig, ax = plt.subplots(figsize=(10, 10))
    divider = make_axes_locatable(ax)
    cbar_axes = [divider.append_axes("bottom", size="3%", pad=0.4 + i * 0.05) for i in range(len(gene_population_pairs))]

    other_cells = subset.obs[~subset.obs['Merged_Cluster_Name_scevan'].isin([pop for _, pop in gene_population_pairs])]
    ax.scatter(
        other_cells['x'], other_cells['y'],
        c="#D3D3D3", s=30, alpha=0.4, edgecolors='none', label='Other cells', zorder=1
    )

    gene_colormaps = {"CD86": "Blues", "CTLA4": "Reds"}

    for i, (gene, population) in enumerate(gene_population_pairs):
        population_data = subset.obs[subset.obs['Merged_Cluster_Name_scevan'] == population][['x', 'y']]
        expression_values = subset[subset.obs['Merged_Cluster_Name_scevan'] == population].to_df()[gene]

        sc = ax.scatter(
            population_data['x'], population_data['y'],
            c=expression_values, cmap=gene_colormaps.get(gene, 'viridis'),
            s=50, alpha=0.9, edgecolors='black', linewidths=0.5,
            label=f"{gene} in {population}", zorder=2
        )

        cbar = plt.colorbar(sc, cax=cbar_axes[i], orientation='horizontal')
        cbar.set_label(f"{gene} Expression", fontsize=10, labelpad=2)
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_aspect(0.2)

    if kde_path is not None:
        contour_patch = mpatches.PathPatch(
            kde_path, facecolor='none', edgecolor='black',
            linewidth=2, linestyle='--', label='KDE-defined region', zorder=3
        )
        ax.add_patch(contour_patch)

    ax.set_title(f"Gene Expression in {sample_id}", fontsize=16, pad=15)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Spatial 1", fontsize=12)
    ax.set_ylabel("Spatial 2", fontsize=12)

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    scale_length = 500
    start_x = x_min + 0.05 * (x_max - x_min)
    start_y = y_min + 0.08 * (y_max - y_min)
    ax.plot([start_x, start_x + scale_length], [start_y, start_y], 'k-', linewidth=2, zorder=4)
    ax.text(start_x + scale_length / 2, start_y - 0.02 * (y_max - y_min),
            f'{scale_length} µm', ha='center', va='top', fontsize=10)

    legend = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9, frameon=True)
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(1.2)

    for spine in ax.spines.values():
        spine.set_linewidth(1.5)

    pdf_filename = f"{save_path}/spatial_gene_expression_{sample_id}.pdf"
    # plt.savefig(pdf_filename, dpi=1200, bbox_inches='tight')
    plt.show()
    plt.close(fig)
    print(f"✅ Saved: {pdf_filename}")


# Example KDE and plot invocation:
# kde_path, x_center, y_center = generate_kde_contour(
#     adata=adata_merged_filtered_seurat,
#     sample_id='your_sample_id',
#     categories=['NSUN5_MUT', 'NSUN5_specific_TCR']
# )
# gene_population_pairs = [
#     ("CD86", "Macrophage"),
#     ("CTLA4", "CD8 T Ex")
# ]
# plot_gene_expression_with_kde(
#     adata=adata_merged_filtered_seurat,
#     sample_id='your_sample_id',
#     gene_population_pairs=gene_population_pairs,
#     kde_path=kde_path
# )
