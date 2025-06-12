# Author: Amit Sud  
# Description: This script calculates Moran’s I for spatial autocorrelation of module scores in tumor cells across samples.
# Input:
#   - AnnData object with tumor cell annotation and spatial coordinates ('x', 'y')
#   - Tab-delimited module file (genes per module, one module per column)
# Output:
#   - DataFrame of Moran’s I values and p-values for each module-sample pair
#   - Top N spatially autocorrelated modules per sample
#   - Optional CSV export of results

import pandas as pd
import numpy as np
import scanpy as sc
from libpysal.weights import KNN
from esda.moran import Moran
from tqdm import tqdm

# --- Parameters ---
module_file = "path/to/ITH_scores.txt"  # Replace with actual file path
k_neighbors = 25
top_n = 20

# --- Step 1: Load module gene sets ---
print(f"Reading module scores from: {module_file}")
module_df = pd.read_csv(module_file, sep='\t')
module_dict = {col: module_df[col].dropna().tolist() for col in module_df.columns}
print(f"Loaded {len(module_dict)} modules.")

# --- Step 2: Initialize results container ---
results = []

# --- Step 3: Iterate over unique samples ---
unique_samples = adata_merged_filtered_seurat.obs["sample_id"].unique()
print(f"Processing {len(unique_samples)} samples: {unique_samples}")

for sample_id in tqdm(unique_samples, desc="Samples"):
    print(f"\nProcessing sample: {sample_id}")
    sample_adata = adata_merged_filtered_seurat[adata_merged_filtered_seurat.obs["sample_id"] == sample_id]

    # Tumor cell subset
    tumor_cells = sample_adata[sample_adata.obs["manual_celltype_annotation"] == "Tumor"].copy()
    print(f"  Tumor cells: {tumor_cells.shape[0]}")

    if tumor_cells.shape[0] < k_neighbors + 1:
        print(f"  Skipping sample {sample_id}: not enough tumor cells for k={k_neighbors}")
        continue

    if not {"x", "y"}.issubset(tumor_cells.obs.columns):
        raise ValueError("Spatial coordinates 'x' and 'y' are missing.")

    coords = tumor_cells.obs[["x", "y"]].values

    # --- Step 4: Compute module scores and Moran's I ---
    for module_name, gene_list in module_dict.items():
        valid_genes = [gene for gene in gene_list if gene in tumor_cells.var_names]
        if len(valid_genes) == 0:
            print(f"  Skipping {module_name}: no valid genes in dataset.")
            continue

        # Score genes
        try:
            sc.tl.score_genes(tumor_cells, gene_list=valid_genes, score_name=module_name)
        except Exception as e:
            print(f"  Error computing score for {module_name}: {e}")
            continue

        values = tumor_cells.obs[module_name].values
        if np.var(values) == 0:
            print(f"  Skipping {module_name}: zero variance.")
            continue

        # Moran's I
        try:
            knn_graph = KNN.from_array(coords, k=k_neighbors)
            moran = Moran(values, knn_graph)
            results.append({
                "sample_id": sample_id,
                "module": module_name,
                "moran_I": moran.I,
                "p_value": moran.p_norm
            })
        except Exception as e:
            print(f"  Moran's I error for {module_name}: {e}")
            continue

# --- Step 5: Create results DataFrame ---
results_df = pd.DataFrame(results)
print(f"\nTotal successful computations: {results_df.shape[0]}")

# --- Step 6: Extract top N per sample ---
top_df = (
    results_df.sort_values(["sample_id", "moran_I"], ascending=[True, False])
    .groupby("sample_id")
    .head(top_n)
    .reset_index(drop=True)
)

print("\nTop results per sample:")
print(top_df)

# --- Optional: Save outputs ---
# results_df.to_csv("path/to/all_morans_results.csv", index=False)
# top_df.to_csv("path/to/top_morans_per_sample.csv", index=False)
