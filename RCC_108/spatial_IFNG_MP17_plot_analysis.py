# Author: Amit Sud
# Description: This script calculates and compares spatial expression patterns of IFNG and a specified
#              tumor module (e.g., interferon response) using kernel density estimation (KDE). It visualizes
#              the spatial overlap and computes Jaccard index, Spearman, and Pearson correlations.
# Input:
#   - AnnData object (adata_merged_filtered_seurat) with:
#       - 'sample_id', 'x', 'y', 'manual_celltype_annotation', 'IFNG_expr', and '[module_name]_score' in obs
#   - Module score file (TSV format) with column 'MP17_Interferon_MHC-II_I'
# Output:
#   - KDE overlap plots and Jaccard/Correlation statistics printed to console

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.stats import gaussian_kde, spearmanr, pearsonr
from skimage.filters import threshold_otsu

def compute_module_score_tumor(adata, module_file, module_name):
    print(f"\n🔍 Reading module score file: {module_file}...\n")
    module_df = pd.read_csv(module_file, sep="\t")

    if module_name not in module_df.columns:
        raise ValueError(f"Module '{module_name}' not found in module score file.")

    module_genes = module_df[module_name].dropna().values
    print(f"✅ Extracted {len(module_genes)} genes for module '{module_name}'.\n")
    
    unique_samples = adata.obs['sample_id'].unique()
    print(f"🔍 Computing module scores independently for {len(unique_samples)} samples...\n")
    
    for sample in unique_samples:
        print(f"➡ Processing sample: {sample}...")
        sample_data = adata[
            (adata.obs['sample_id'] == sample) &
            (adata.obs['manual_celltype_annotation'] == 'Tumor')
        ].copy()
        
        if sample_data.n_obs == 0:
            print(f"⚠ No Tumor cells found in sample {sample}. Skipping.")
            continue

        sc.tl.score_genes(sample_data, gene_list=module_genes, score_name=f"{module_name}_score")
        adata.obs.loc[sample_data.obs.index, f"{module_name}_score"] = sample_data.obs[f"{module_name}_score"]
        print(f"✅ Module score added for sample {sample}.")
    
    return adata

def compute_kde(adata, sample_id, feature_column, tumor_only=False, grid_size=250):
    if tumor_only:
        sample_data = adata[
            (adata.obs['sample_id'] == sample_id) &
            (adata.obs['manual_celltype_annotation'] == 'Tumor')
        ]
    else:
        sample_data = adata[adata.obs['sample_id'] == sample_id]
    
    x = sample_data.obs["x"].values
    y = sample_data.obs["y"].values
    values = sample_data.obs[feature_column].values
    
    nonzero_mask = values > 0
    x_filtered, y_filtered, values_filtered = x[nonzero_mask], y[nonzero_mask], values[nonzero_mask]
    
    if len(x_filtered) == 0:
        raise ValueError(f"No nonzero values found for {feature_column} in sample {sample_id}.")
    
    n, d = len(x_filtered), 2
    bandwidth = (n * (d + 2) / 4) ** (-1 / (d + 4))
    
    xy = np.vstack([x_filtered, y_filtered])
    kde = gaussian_kde(xy, weights=values_filtered, bw_method=bandwidth)
    
    x_min, x_max = x_filtered.min(), x_filtered.max()
    y_min, y_max = y_filtered.min(), y_filtered.max()
    x_grid, y_grid = np.meshgrid(
        np.linspace(x_min, x_max, grid_size),
        np.linspace(y_min, y_max, grid_size)
    )
    kde_density = kde(np.vstack([x_grid.ravel(), y_grid.ravel()])).reshape(x_grid.shape)
    
    return x_grid, y_grid, kde_density

# Example usage (replace with appropriate file and module name)
module_file_path = "path/to/module_score_file.tsv"
module_name = "MP17_Interferon_MHC-II_I"

adata_merged_filtered_seurat = compute_module_score_tumor(
    adata_merged_filtered_seurat,
    module_file=module_file_path,
    module_name=module_name
)

# Compute KDEs
x_grid, y_grid, kde_ifng = compute_kde(adata_merged_filtered_seurat, 'rcc_108_2', 'IFNG_expr', tumor_only=False)
_, _, kde_mp17 = compute_kde(adata_merged_filtered_seurat, 'rcc_108_2', f'{module_name}_score', tumor_only=True)

# Plot KDE Overlap
fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].imshow(kde_ifng, origin='lower', cmap='Blues', extent=[x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()])
ax[0].set_title("IFNG KDE Density (All Cells)")
ax[1].imshow(kde_mp17, origin='lower', cmap='Reds', extent=[x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()])
ax[1].set_title(f"{module_name} KDE Density (Tumor Cells)")

combined_kde = (kde_ifng + kde_mp17) / 2
ax[2].imshow(combined_kde, origin='lower', cmap='Purples', extent=[x_grid.min(), x_grid.max(), y_grid.min(), y_grid.max()])
ax[2].set_title("Overlap of KDE Densities")

plt.tight_layout()
plt.show()

# Jaccard Index
threshold_ifng = threshold_otsu(kde_ifng)
threshold_mp17 = threshold_otsu(kde_mp17)
binary_ifng = kde_ifng > threshold_ifng
binary_mp17 = kde_mp17 > threshold_mp17
intersection = np.logical_and(binary_ifng, binary_mp17).sum()
union = np.logical_or(binary_ifng, binary_mp17).sum()
jaccard_index = intersection / union
print(f"🔍 Jaccard Index (Hotspot Overlap): {jaccard_index:.3f}")

# Spearman Correlation
spearman_corr, spearman_pval = spearmanr(kde_ifng.flatten(), kde_mp17.flatten())
print(f"🔍 Spearman Correlation: {spearman_corr:.3f} (p={spearman_pval:.3e})")

# Pearson Correlation (log-transformed)
log_kde_ifng = np.log1p(kde_ifng.flatten())
log_kde_mp17 = np.log1p(kde_mp17.flatten())
pearson_corr, pearson_pval = pearsonr(log_kde_ifng, log_kde_mp17)
print(f"🔍 Pearson Correlation (log-transformed): {pearson_corr:.3f} (p={pearson_pval:.3e})")
