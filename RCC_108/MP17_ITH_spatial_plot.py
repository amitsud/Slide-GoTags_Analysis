# Author: Amit Sud
# Description: This script computes module scores for tumor cells across samples,
#              and plots a kernel density estimate (KDE) of spatial distribution
#              of the scores, clipped by a user-defined polygon.

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from shapely.geometry import Point, Polygon

def compute_module_score_tumor(adata, module_file, module_name):
    """
    Computes the module score for Tumor cells for each unique sample.

    Parameters:
    - adata: AnnData object
    - module_file: Path to tab-delimited file with gene modules
    - module_name: Column name in module_file representing the desired gene module

    Returns:
    - Updated AnnData object with module scores added to .obs
    """
    print(f"\n🔍 Reading module score file: {module_file}...\n")
    module_df = pd.read_csv(module_file, sep="\t")

    if module_name not in module_df.columns:
        raise ValueError(f"Module '{module_name}' not found in module score file.")

    module_genes = module_df[module_name].dropna().values
    print(f"✅ Extracted {len(module_genes)} genes for module '{module_name}'.\n")
    
    score_column = f"{module_name}_score"
    adata.obs[score_column] = np.nan

    unique_samples = adata.obs['sample_id'].unique()
    print(f"🔍 Computing module scores independently for {len(unique_samples)} samples...\n")

    for sample in unique_samples:
        print(f"➡ Processing sample: {sample}...")
        sample_data = adata[(adata.obs['sample_id'] == sample) & (adata.obs['manual_celltype_annotation'] == 'Tumor')].copy()

        if sample_data.n_obs == 0:
            print(f"⚠ No Tumor cells found in sample {sample}. Skipping.")
            continue

        sc.tl.score_genes(sample_data, gene_list=module_genes, score_name=score_column)
        adata.obs.loc[sample_data.obs.index, score_column] = sample_data.obs[score_column]
        print(f"✅ Module score added for sample {sample}.")

    return adata


def plot_kde_with_polygon_clipping(adata, sample_id, module_name, polygon_coords, grid_size=250, percentile=90, save_pdf=False, pdf_filename="kde_plot_with_polygon.pdf"):
    """
    Plots a KDE for the module score with polygon clipping.

    Parameters:
    - adata: AnnData object
    - sample_id: Sample to plot
    - module_name: Module score column name (without _score suffix)
    - polygon_coords: Polygon coordinates for clipping
    - grid_size: Resolution of the KDE grid
    - percentile: Threshold percentile for contour
    - save_pdf: Whether to save the plot
    - pdf_filename: Output filename for the plot
    """
    sample_data = adata[adata.obs['sample_id'] == sample_id].copy()
    score_column = f"{module_name}_score"

    if score_column not in sample_data.obs.columns:
        raise ValueError(f"❌ ERROR: Module score '{score_column}' not found in .obs!")

    valid_mask = sample_data.obs[score_column] > 0
    x = sample_data.obs.loc[valid_mask, 'x']
    y = sample_data.obs.loc[valid_mask, 'y']
    scores = sample_data.obs.loc[valid_mask, score_column]

    if len(x) == 0:
        print(f"⚠ No nonzero module scores found in sample {sample_id}.")
        return

    kde = gaussian_kde(np.vstack([x, y]), weights=scores, bw_method='scott')

    x_min, x_max, y_min, y_max = x.min(), x.max(), y.min(), y.max()
    x_grid, y_grid = np.meshgrid(np.linspace(x_min, x_max, grid_size),
                                 np.linspace(y_min, y_max, grid_size))
    grid_coords = np.vstack([x_grid.ravel(), y_grid.ravel()])
    density = kde(grid_coords).reshape(grid_size, grid_size)

    polygon = Polygon(polygon_coords)
    mask = np.array([[polygon.contains(Point(xi, yi)) for xi, yi in zip(row_x, row_y)] for row_x, row_y in zip(x_grid, y_grid)])
    density[~mask] = np.nan

    threshold = np.nanpercentile(density, percentile)

    fig, ax = plt.subplots(figsize=(10, 8))
    kde_plot = ax.imshow(density, origin='lower',
                         extent=(x_min, x_max, y_min, y_max),
                         cmap='coolwarm', alpha=0.8)

    poly_x, poly_y = zip(*polygon_coords)
    ax.plot(poly_x + (poly_x[0],), poly_y + (poly_y[0],), 'white', linewidth=2)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)

    bottom_right_x = x_min + 0.90 * (x_max - x_min)
    bottom_right_y = y_min + 0.08 * (y_max - y_min)
    size_microns = 500
    ax.plot([bottom_right_x, bottom_right_x - size_microns], [bottom_right_y, bottom_right_y], 'k-', linewidth=2)
    ax.text(bottom_right_x - size_microns / 2, bottom_right_y - 0.015 * (y_max - y_min),
            f'{size_microns} µm', ha='center', fontsize=12, color='black')

    plt.title(f"KDE of {module_name} Score (Sample {sample_id})", fontsize=14)
    plt.colorbar(kde_plot, ax=ax, label=f"{module_name} Score Density")

    if save_pdf:
        plt.savefig(pdf_filename, dpi=1200, bbox_inches='tight', format='pdf')
        print(f"✅ Plot saved as: {pdf_filename}")

    plt.show()


# ✅ Example Usage

# Step 1: Compute module scores
adata_merged_filtered_seurat = compute_module_score_tumor(
    adata=adata_merged_filtered_seurat,
    module_file='path/to/module_scores.txt',
    module_name='MP17_Interferon_MHC-II_I'
)

# Step 2: Plot KDE with polygon clipping
plot_kde_with_polygon_clipping(
    adata=adata_merged_filtered_seurat,
    sample_id='rcc_108_2',
    module_name='MP17_Interferon_MHC-II_I',
    polygon_coords=polygon_coords,
    grid_size=250,
    percentile=90,
    save_pdf=True,
    pdf_filename="MP17_kde_plot_with_clipping.pdf"
)

