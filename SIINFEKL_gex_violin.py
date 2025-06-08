# Author: Amit Sud
# Date: 1st May 2025
# Description: Plots a violin plot showing the normalized expression of the gene SIINFEKL_WPRE across annotated cell types.
#              It uses raw UMI counts to compute the median expression where the gene is detected, and plots normalized values.
# Input:
#   - adata_filtered: AnnData object with:
#       • adata_filtered.var_names including 'SIINFEKL_WPRE'
#       • adata_filtered.layers['counts'] for raw UMI counts
#       • adata_filtered.X for normalized expression
#       • adata_filtered.obs['manual_celltype_annotation'] for cell type labels
# Output:
#   - A violin plot saved as "SIINFEKL_WPRE_violin_plot.pdf"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse

# Define a rainbow color palette using seaborn and custom hex codes
rainbow_palette = sns.color_palette("hsv", 7).as_hex()
color_palette = {
    "T-cells": '#FF0000',              # Red
    "Tumor cells": '#333399',          # Blue
    "Myeloid cells": '#9B30FF',        # Bright Purple
    "Fibroblasts": rainbow_palette[3]  # Greenish tone from rainbow palette
}

# Ensure the gene of interest is present in the dataset
if 'SIINFEKL_WPRE' in adata_filtered.var_names:
    # Check for raw count data
    if 'counts' in adata_filtered.layers:
        gene_expression_raw = adata_filtered[:, 'SIINFEKL_WPRE'].layers['counts']
        
        # Convert to dense array if sparse
        if sparse.issparse(gene_expression_raw):
            gene_expression_raw = gene_expression_raw.toarray().flatten()
        else:
            gene_expression_raw = gene_expression_raw.flatten()
        
        # Filter out cells with 0 UMI for the gene
        detected_umi_counts = gene_expression_raw[gene_expression_raw > 0]
        median_umi_detected = np.median(detected_umi_counts)
        print(f"The median number of UMI detected across cells for SIINFEKL_WPRE (where detected) is: {median_umi_detected}")
    else:
        print("Counts layer not available in the AnnData object.")
    
    # Extract normalized expression (e.g., log-normalized or scaled)
    gene_expression = adata_filtered[:, 'SIINFEKL_WPRE'].X
    if sparse.issparse(gene_expression):
        gene_expression = gene_expression.toarray().flatten()
    else:
        gene_expression = gene_expression.flatten()

    # Extract cell type labels
    cell_types = adata_filtered.obs['manual_celltype_annotation']

    # Build DataFrame for plotting
    data = pd.DataFrame({
        'SIINFEKL_WPRE': gene_expression,
        'Cell Type': cell_types
    })

    # Generate custom palette for the available cell types
    custom_palette = {ct: color_palette[ct] for ct in data['Cell Type'].unique() if ct in color_palette}

    # Plot violin plot
    plt.figure(figsize=(12, 8))
    sns.violinplot(
        x='Cell Type',
        y='SIINFEKL_WPRE',
        data=data,
        inner='quartile',
        palette=custom_palette
    )
    plt.title('Violin plot of normalized gene expression of SIINFEKL_WPRE by cell type', fontsize=16)
    plt.xlabel('Cell Type', fontsize=14)
    plt.ylabel('Normalized Expression of SIINFEKL_WPRE', fontsize=14)
    plt.xticks(rotation=45, fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(False)

    # Save the figure
    plt.savefig("SIINFEKL_WPRE_violin_plot.pdf", format="pdf", bbox_inches="tight", dpi=1200)
    plt.show()

else:
    print("Gene 'SIINFEKL_WPRE' not found in the dataset.")
