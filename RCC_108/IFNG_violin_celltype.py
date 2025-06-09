# Author: Amit Sud
# Date: 1st May 2025
# Description: This script visualizes IFNG gene expression across annotated cell types 
#              in a specific RCC sample using a violin plot.
# Input:
#   - AnnData object (adata_merged_filtered_seurat) with:
#       - 'sample_id' == 'rcc_108_2'
#       - 'manual_celltype_annotation' in obs
#       - 'IFNG' in var_names
# Output:
#   - Violin plot of IFNG expression by cell type, displayed via matplotlib

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# --- Parameters ---
sample_id = 'rcc_108_2'
gene = 'IFNG'

# --- Color palette ---
color_palette = {
    "Tumor": "#A5CEE3",
    "NK cell": "#F7DC6F",
    "Fibroblast": "#B1CC70",
    "T-cell": "#E32128",
    "Endothelial cell": "#F38DB2",
    "Myeloid": "#F57F20",
    "B-cell": "#6A0DAD",
}

# --- Step 1: Filter for sample ---
adata_filtered = adata_merged_filtered_seurat[
    adata_merged_filtered_seurat.obs['sample_id'] == sample_id
].copy()

# --- Step 2: Check for gene ---
if gene not in adata_filtered.var_names:
    raise ValueError(f"Gene {gene} not found in the dataset.")

# --- Step 3: Extract expression values ---
expr = adata_filtered[:, gene].X
expr = expr.toarray().flatten() if not isinstance(expr, np.ndarray) else expr.flatten()

# --- Step 4: Prepare plotting data ---
df = pd.DataFrame({
    f'{gene}_expression': expr,
    'manual_celltype_annotation': adata_filtered.obs['manual_celltype_annotation'].values
})

# --- Step 5: Set color palette based on observed cell types ---
unique_cell_types = df['manual_celltype_annotation'].unique()
palette = {ct: color_palette.get(ct, "gray") for ct in unique_cell_types}

# --- Step 6: Plot ---
plt.figure(figsize=(10, 6))
sns.violinplot(
    data=df,
    x='manual_celltype_annotation',
    y=f'{gene}_expression',
    scale='width',
    palette=palette
)
plt.xticks(rotation=45, ha='right')
plt.xlabel('Manual Cell Type Annotation')
plt.ylabel(f'{gene} Expression')
plt.title(f'{gene} Expression by Cell Type in {sample_id}')
plt.tight_layout()
plt.show()

