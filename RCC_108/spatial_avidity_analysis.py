# Author: Amit Sud
# Date: 8th June 2025
# Description: This script maps NSUN5-specific TCR avidity scores (CDR3A and CDR3B) to cells in a spatial
#              single-cell AnnData object and tests whether low-avidity TCRs are spatially enriched in
#              a defined region (via Fisher's Exact Test).
# Input:
#   - TCR avidity table (tab-delimited) with columns: 'CDR3A', 'CDR3B', 'NSUN5_Avidity', 'Specificity'
#   - AnnData object (adata_merged_filtered_seurat) with:
#       - 'sample_id', 'inside_region', 'TRA_cdr3', 'TRB_cdr3' in .obs
# Output:
#   - Updated AnnData object with columns:
#       - 'NSUN5_avidity_TRA_cdr3' and 'NSUN5_avidity_TRB_cdr3'
#   - Printed summary of matched CDR3s and results of Fisher's Exact Test for spatial enrichment

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

# ✅ Step 1: Read NSUN5-specific TCR avidity file
file_path = "path/to/108_TCR_avidity.txt"
print(f"🔍 Reading file: {file_path}")
rcc_tcr_specificity = pd.read_csv(file_path, sep="\t")

# ✅ Step 2: Validate expected columns
expected_columns = ["CDR3A", "CDR3B", "NSUN5_Avidity", "Specificity"]
for col in expected_columns:
    assert col in rcc_tcr_specificity.columns, f"⚠ Column '{col}' missing in avidity data"
print("✅ All required columns are present.")

# ✅ Step 3: Filter for NSUN5-specific TCRs
rcc_tcr_specificity_filtered = rcc_tcr_specificity[
    rcc_tcr_specificity["Specificity"] == "NSUN5 p.Q18K"
]
print(f"✅ Filtered for Specificity = 'NSUN5 p.Q18K': {rcc_tcr_specificity_filtered.shape[0]} rows retained.")

# ✅ Step 4: Create mapping dictionaries
cdr3a_to_avidity = dict(zip(rcc_tcr_specificity_filtered["CDR3A"], rcc_tcr_specificity_filtered["NSUN5_Avidity"]))
cdr3b_to_avidity = dict(zip(rcc_tcr_specificity_filtered["CDR3B"], rcc_tcr_specificity_filtered["NSUN5_Avidity"]))
print(f"✅ Created mapping dictionaries: {len(cdr3a_to_avidity)} TRA, {len(cdr3b_to_avidity)} TRB")

# ✅ Step 5: Map avidity values to AnnData object
adata_merged_filtered_seurat.obs["NSUN5_avidity_TRA_cdr3"] = adata_merged_filtered_seurat.obs["TRA_cdr3"].map(cdr3a_to_avidity)
adata_merged_filtered_seurat.obs["NSUN5_avidity_TRB_cdr3"] = adata_merged_filtered_seurat.obs["TRB_cdr3"].map(cdr3b_to_avidity)

# ✅ Step 6: Summary statistics
print("\n🔍 Summary of Mapped Avidity Scores:")
print(f"✔ TRA matches: {adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].notna().sum()}")
print(f"✔ TRB matches: {adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].notna().sum()}")

# ✅ Step 7: Random sample preview
print("\n🔍 Sample TRA mappings:")
print(adata_merged_filtered_seurat.obs[['TRA_cdr3', 'NSUN5_avidity_TRA_cdr3']].dropna().sample(5))

print("\n🔍 Sample TRB mappings:")
print(adata_merged_filtered_seurat.obs[['TRB_cdr3', 'NSUN5_avidity_TRB_cdr3']].dropna().sample(5))

# ✅ Step 8: Filter for sample of interest
sample_id = 'rcc_108_2'
rcc_108_2_data = adata_merged_filtered_seurat[adata_merged_filtered_seurat.obs['sample_id'] == sample_id].obs.copy()

# ✅ Step 9: Define Fisher’s Exact Test function
def perform_fishers_exact(data, column, threshold=80):
    """
    Perform Fisher's Exact Test to assess spatial enrichment of low-avidity TCRs.

    Parameters:
    - data: DataFrame filtered for a specific sample.
    - column: Name of the column containing avidity scores.
    - threshold: Avidity threshold to define 'low' vs 'high' (default = 80).

    Returns:
    - odds_ratio, p_value
    """
    data_filtered = data.dropna(subset=[column])

    inside_true = (data_filtered['inside_region'] == True) & (data_filtered[column] < threshold)
    inside_false = (data_filtered['inside_region'] == True) & (data_filtered[column] >= threshold)
    outside_true = (data_filtered['inside_region'] == False) & (data_filtered[column] < threshold)
    outside_false = (data_filtered['inside_region'] == False) & (data_filtered[column] >= threshold)

    contingency_table = [
        [inside_true.sum(), inside_false.sum()],
        [outside_true.sum(), outside_false.sum()]
    ]

    odds_ratio, p_value = fisher_exact(contingency_table)

    print(f"\n🔍 Fisher's Exact Test for {column}:")
    print(f"📊 Contingency Table: {contingency_table}")
    print(f"✔ Odds Ratio: {odds_ratio:.3f}")
    print(f"✔ P-value: {p_value:.2e}")

    return odds_ratio, p_value

# ✅ Step 10: Run Fisher’s test for TRA and TRB
odds_ratio_TRA, p_value_TRA = perform_fishers_exact(rcc_108_2_data, 'NSUN5_avidity_TRA_cdr3')
odds_ratio_TRB, p_value_TRB = perform_fishers_exact(rcc_108_2_data, 'NSUN5_avidity_TRB_cdr3')
