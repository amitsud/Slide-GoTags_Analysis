# Author: Amit Sud  
# Date: 1st May 2025  
# Description: This script calculates module scores in T-cells, computes their distances to Rpl18_MUT tumor cells,
#              stratifies cells by clonotype and proximity, and visualizes results using a heatmap and scatter plots.
# Input:
#   - AnnData object: adata_merged_seurat_filtered
#   - Module score file: tab-delimited text file of gene modules
#   - Required columns in .obs:
#       • 'sample_id', 'manual_celltype_annotation', 'Rpl18_mut_tcr', 'Rpl18_TCR', 'x', 'y'
# Output:
#   - Heatmap PDF and per-module scatter plot PDFs showing distance vs. module score relationships

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import linregress
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap

# === Step 1: Load Module Scores ===
module_path = "/n/data2/dfci/medonc/cwu/amit/reference/module_scores/Mouse module_T_PD1 1-14-2025.txt"
module_scores = pd.read_csv(module_path, delimiter='\t').dropna(how='all')

# === Step 2: Subset adata ===
sample_id = 'mc38_t1_pd1'
filtered_adata = adata_merged_seurat_filtered[
    adata_merged_seurat_filtered.obs['sample_id'] == sample_id
]

if filtered_adata.shape[0] == 0:
    raise ValueError(f"No data found for sample {sample_id}.")

# === Step 3: Subset to T-cells ===
t_cell_data = filtered_adata[filtered_adata.obs['manual_celltype_annotation'] == 'T-cell']

if t_cell_data.shape[0] == 0:
    raise ValueError(f"No T-cells found in sample {sample_id}.")

# === Step 4: Get coordinates of Rpl18_MUT tumor cells ===
mut_coords = filtered_adata.obs[
    filtered_adata.obs['Rpl18_mut_tcr'] == 'Rpl18_MUT'
][['x', 'y']].to_numpy()

if len(mut_coords) == 0:
    raise ValueError(f"No Rpl18_MUT cells found for distance calculations.")

# === Step 5: Calculate Module Scores ===
for module_name in module_scores.columns:
    gene_list = module_scores[module_name].dropna().tolist()
    if gene_list:
        sc.tl.score_genes(t_cell_data, gene_list=gene_list, score_name=f"{module_name}_score")

module_score_columns = [f"{col}_score" for col in module_scores.columns]
print("Module scores calculated:", module_score_columns)

# === Step 6: Calculate Nearest Distance to Rpl18_MUT ===
def calculate_nearest_distances(source_coords, target_coords):
    if len(source_coords) == 0 or len(target_coords) == 0:
        return np.array([])
    tree = cKDTree(target_coords)
    distances, _ = tree.query(source_coords, k=1)
    return distances

specific_coords = t_cell_data.obs[t_cell_data.obs['Rpl18_TCR'] == 'Rpl18_specific_TCR'][['x', 'y']].to_numpy()
other_coords = t_cell_data.obs[t_cell_data.obs['Rpl18_TCR'] == 'Other_clonotype'][['x', 'y']].to_numpy()

specific_dists = calculate_nearest_distances(specific_coords, mut_coords)
other_dists = calculate_nearest_distances(other_coords, mut_coords)

t_cell_data.obs['Distance'] = np.nan
if len(specific_dists):
    t_cell_data.obs.loc[t_cell_data.obs['Rpl18_TCR'] == 'Rpl18_specific_TCR', 'Distance'] = specific_dists
if len(other_dists):
    t_cell_data.obs.loc[t_cell_data.obs['Rpl18_TCR'] == 'Other_clonotype', 'Distance'] = other_dists

t_cell_data.obs['Distance_Group'] = np.where(t_cell_data.obs['Distance'] < 50, '<50', '>=50')

# === Step 7: Prepare Heatmap Data ===
heatmap_data = []
cell_counts = {}

for clonotype in ['Rpl18_specific_TCR', 'Other_clonotype']:
    for group in ['<50', '>=50']:
        subset = t_cell_data.obs[
            (t_cell_data.obs['Rpl18_TCR'] == clonotype) &
            (t_cell_data.obs['Distance_Group'] == group)
        ]
        if not subset.empty:
            avg_scores = subset[module_score_columns].mean().to_frame(name=f"{clonotype}_{group}")
            heatmap_data.append(avg_scores)
            cell_counts[f"{clonotype}_{group}"] = len(subset)

if not heatmap_data:
    raise ValueError("No valid data for heatmap.")
heatmap_df = pd.concat(heatmap_data, axis=1)

# === Step 8: Plot Heatmap ===
cmap = LinearSegmentedColormap.from_list("custom", ['#6eadd9', 'white', '#df88a4'])

plt.figure(figsize=(12, 12))
sns.heatmap(
    heatmap_df, annot=True, fmt=".2f", cmap=cmap, cbar_kws={'label': 'Mean module score'},
    linewidths=0.8, linecolor='black', square=True, annot_kws={"size": 9}
)

legend_text = "\n".join([f"{k}: {v} cells" for k, v in cell_counts.items()])
plt.gcf().text(1.02, 0.5, legend_text, fontsize=10, va='center')
plt.title("Heatmap of Mean Module Scores (T-cells stratified by distance to Rpl18_MUT)")
plt.ylabel("Modules")
plt.xlabel("Clonotype + Distance Group")
plt.xticks(rotation=45)
plt.tight_layout()

plt.savefig("t_cell_module_scores_heatmap.pdf", format="pdf", dpi=1200)
plt.show()

# === Step 9: Scatter Plots with Linear Regression ===
for module_name in module_scores.columns:
    module_col = f"{module_name}_score"
    plt.figure(figsize=(8, 6))

    for clonotype, color in [('Rpl18_specific_TCR', '#FF0000'), ('Other_clonotype', '#FFCCCC')]:
        subset = t_cell_data.obs[t_cell_data.obs['Rpl18_TCR'] == clonotype]
        x = subset['Distance'].to_numpy()
        y = subset[module_col].to_numpy()

        if len(x) == 0 or len(y) == 0:
            print(f"Skipping {clonotype}, Module: {module_name} due to insufficient data.")
            continue

        sorted_idx = np.argsort(x)
        x_sorted, y_sorted = x[sorted_idx], y[sorted_idx]

        slope, intercept, r, p, stderr = linregress(x_sorted, y_sorted)

        plt.scatter(x_sorted, y_sorted, s=100, alpha=0.9, color=color, label=clonotype, edgecolor='black', linewidth=0.5)
        plt.plot(x_sorted, intercept + slope * x_sorted, color=color, label=f"{clonotype} (p = {p:.4f})")

    plt.xlabel("Distance to Rpl18_MUT")
    plt.ylabel(f"{module_name} Score")
    plt.title(f"{module_name} Score vs. Distance to Rpl18_MUT")
    plt.legend(fontsize=10)
    plt.grid(False)
    plt.tight_layout()

    plt.savefig(f"scatter_{module_name}_distance.pdf", format="pdf", dpi=1200)
    plt.show()
