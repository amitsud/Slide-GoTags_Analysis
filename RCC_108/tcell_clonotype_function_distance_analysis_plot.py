# nearest_distance_module_score_analysis.py
# Author: Your Name
# Description: Computes nearest distances between CD8 T cells and target cells (e.g., NSUN5_MUT),
# assigns them to categories, calculates module scores, and performs linear regression
# between distance and score for various clonotypes.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import linregress
import scanpy as sc

# Placeholder: Replace with your loaded AnnData object and module scores
adata = YOUR_ANN_DATA_OBJECT
module_scores = YOUR_MODULE_SCORE_TABLE  # A DataFrame with columns of module names and gene lists

# Lists to store results
all_combined_data_nsun5, all_combined_data_tumor, all_combined_data_other = [], [], []
significant_modules = []

for sample in adata.obs['sample_id'].unique():
    print(f"Processing sample: {sample}")
    sample_data = adata[adata.obs['sample_id'] == sample]

    cd8_nsun5_data = sample_data[
        (sample_data.obs['celltype'].str.contains('CD8 T', case=False, na=False)) &
        (sample_data.obs['TCR_class'] == 'NSUN5_specific_TCR')
    ]
    cd8_tumor_data = sample_data[
        (sample_data.obs['celltype'].str.contains('CD8 T', case=False, na=False)) &
        (sample_data.obs['TCR_class'] == 'Tumor_specific_TCR')
    ]
    cd8_other_data = sample_data[
        (sample_data.obs['celltype'].str.contains('CD8 T', case=False, na=False)) &
        (sample_data.obs['TCR_class'] == 'Other_clonotype_TCR')
    ]
    target_cells = sample_data[sample_data.obs['TCR_class'] == 'Target_Cell_Label']

    if target_cells.shape[0] == 0:
        print(f"Skipping sample {sample} (no target cells).")
        continue

    for group_data, label, storage_list in zip(
        [cd8_nsun5_data, cd8_tumor_data, cd8_other_data],
        ['NSUN5_specific_TCR', 'Tumor_specific_TCR', 'Other_clonotype_TCR'],
        [all_combined_data_nsun5, all_combined_data_tumor, all_combined_data_other]
    ):
        if group_data.shape[0] > 0:
            tree = cKDTree(target_cells.obs[['x', 'y']].values)
            distances, _ = tree.query(group_data.obs[['x', 'y']].values)
            group_data.obs['nearest_distance'] = distances
            storage_list.append(group_data.obs[['nearest_distance', 'sample_id']])

# Combine data
combined_data_nsun5 = pd.concat(all_combined_data_nsun5)
combined_data_tumor = pd.concat(all_combined_data_tumor)
combined_data_other = pd.concat(all_combined_data_other)

# Store nearest distances in the original AnnData object
adata.obs['nearest_distance_nsun5'] = combined_data_nsun5['nearest_distance']
adata.obs['nearest_distance_tumor'] = combined_data_tumor['nearest_distance']
adata.obs['nearest_distance_other'] = combined_data_other['nearest_distance']

# Compute module scores and regressions
for module in module_scores.columns:
    print(f"Processing module: {module}")
    genes = module_scores[module].dropna().tolist()
    available_genes = [gene for gene in genes if gene in adata.var_names]

    if not available_genes:
        print(f"No valid genes for module {module}. Skipping.")
        continue

    sc.tl.score_genes(adata, gene_list=available_genes, score_name=f"{module}_score")

    for group_data, label, color, combined_df in zip(
        [combined_data_nsun5, combined_data_tumor, combined_data_other],
        ['NSUN5_specific_TCR', 'Tumor_specific_TCR', 'Other_clonotype_TCR'],
        ['#FF0000', '#8c78c3', '#FFCCCC'],
        [combined_data_nsun5, combined_data_tumor, combined_data_other]
    ):
        combined_df[f"{module}_score"] = adata.obs.loc[combined_df.index, f"{module}_score"]

    # Linear regression
    regressions = []
    for data in [combined_data_nsun5, combined_data_tumor, combined_data_other]:
        x = data['nearest_distance']
        y = data[f"{module}_score"]
        if not x.empty and not y.empty:
            slope, intercept, r, p, _ = linregress(x, y)
        else:
            p = np.nan
        regressions.append((x, y, slope, intercept, p))

    if any(p < 0.05 for *_, p in regressions):
        significant_modules.append(module)

        # Plotting
        plt.figure(figsize=(8, 6))
        for (x, y, slope, intercept, p), label, color in zip(regressions,
                                                             ['NSUN5_specific_TCR', 'Tumor_specific_TCR', 'Other_clonotype_TCR'],
                                                             ['#FF0000', '#8c78c3', '#FFCCCC']):
            if not x.empty:
                x_np, y_np = x.to_numpy(), y.to_numpy()
                plt.scatter(x_np, y_np, alpha=0.7, color=color, edgecolor='black', linewidth=0.5, label=f'{label}')
                plt.plot(x_np, slope * x_np + intercept, color=color, label=f'{label} (p={p:.2g})')

        plt.xlabel("Nearest Distance to Target Cells")
        plt.ylabel(f"{module} Module Score")
        plt.title(f"{module} Score vs. Distance")
        plt.legend()
        plt.tight_layout()
        plt.grid(False)
        plt.show()
