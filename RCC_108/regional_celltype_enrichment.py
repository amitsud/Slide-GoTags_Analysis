# Author: Amit Sud
# Description: This script tests for enrichment of annotated cell types in defined spatial regions
#              using either Fisher’s exact test or the Chi-squared test based on minimum count thresholds.
#              It calculates odds ratios, log2 fold changes, and p-values for each cell type.
# Input:
#   - AnnData object (adata) with obs containing:
#       - 'sample_id', 'manual_celltype_annotation', 'inside_region'
# Output:
#   - DataFrame summarizing enrichment results per cell type:
#       - Odds Ratio, log2 Fold Change, p-value

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency

def test_celltype_enrichment(adata, sample_id=None, method="fisher", min_count=5):
    """
    Test for spatial enrichment of each cell type inside vs. outside a region.

    Parameters:
        adata (AnnData): The AnnData object with cell annotations.
        sample_id (str, optional): Filter by specific sample_id.
        method (str): "fisher" or "chi2" to select statistical test.
        min_count (int): Threshold for using Chi-squared instead of Fisher’s test.

    Returns:
        pd.DataFrame: Summary with odds ratios, log2FC, and p-values for each cell type.
    """
    df = adata.obs.copy()
    if sample_id:
        df = df[df['sample_id'] == sample_id]

    if 'manual_celltype_annotation' not in df.columns or 'inside_region' not in df.columns:
        raise ValueError("Columns 'manual_celltype_annotation' and 'inside_region' must be present in adata.obs.")

    results = []
    for cell_type in df['manual_celltype_annotation'].unique():
        inside_target = ((df['manual_celltype_annotation'] == cell_type) & (df['inside_region'] == True)).sum()
        outside_target = ((df['manual_celltype_annotation'] == cell_type) & (df['inside_region'] == False)).sum()
        inside_other = ((df['manual_celltype_annotation'] != cell_type) & (df['inside_region'] == True)).sum()
        outside_other = ((df['manual_celltype_annotation'] != cell_type) & (df['inside_region'] == False)).sum()

        contingency_table = [[inside_target, outside_target], [inside_other, outside_other]]

        if method == "chi2" and min(inside_target, outside_target, inside_other, outside_other) >= min_count:
            chi2, p_value, _, _ = chi2_contingency(contingency_table)
            odds_ratio = (inside_target / (inside_target + inside_other)) / (outside_target / (outside_target + outside_other))
        else:
            odds_ratio, p_value = fisher_exact(contingency_table)

        inside_freq = inside_target / (inside_target + inside_other) if (inside_target + inside_other) > 0 else 0
        outside_freq = outside_target / (outside_target + outside_other) if (outside_target + outside_other) > 0 else 1e-9
        log2fc = np.log2(inside_freq / outside_freq) if inside_freq > 0 and outside_freq > 0 else np.nan

        results.append({
            'Cell_Type': cell_type,
            'Inside_Count': inside_target,
            'Outside_Count': outside_target,
            'Odds_Ratio': odds_ratio,
            'log2FC': log2fc,
            'P_Value': p_value
        })

    results_df = pd.DataFrame(results).sort_values('P_Value')
    return results_df

# Example usage:
# results_df = test_celltype_enrichment(
#     adata=adata_merged_filtered_seurat,
#     sample_id='your_sample_id',
#     method="fisher"  # or "chi2"
# )
# print(results_df)
