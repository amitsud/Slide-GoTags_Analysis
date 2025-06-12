# Author: Amit Sud
# Description: This script maps NSUN5-specific TCR avidity scores (CDR3A and CDR3B) to cells in a spatial
#              single-cell AnnData object and stratifies T cells into strong/weak categories based on avidity.
#              It then tests whether low-avidity (strong) TCRs are enriched within a spatial region of interest.
# Input:
#   - TCR avidity file (tab-delimited text with columns: CDR3A, CDR3B, NSUN5_Avidity, Specificity)
#   - AnnData object (adata_merged_filtered_seurat) with:
#       - TRA_cdr3, TRB_cdr3, NSUN5_MUT_TCR, inside_region, sample_id, manual_celltype_annotation in obs
# Output:
#   - Updates adata.obs with avidity scores
#   - Printed summary of cell counts, enrichment statistics, and Fisher’s Exact Test results

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

# Step 1: Read in the avidity file
file_path = "path/to/108_TCR_avidity.txt"
print(f"🔍 Reading file: {file_path}")
rcc_tcr_specificity = pd.read_csv(file_path, sep="\t")

# Step 2: Validate columns
expected_columns = ["CDR3A", "CDR3B", "NSUN5_Avidity", "Specificity"]
for col in expected_columns:
    assert col in rcc_tcr_specificity.columns, f"⚠ Column '{col}' missing in avidity data"
print("✅ All required columns are present.")

# Step 3: Filter for NSUN5-specific TCRs
rcc_tcr_specificity_filtered = rcc_tcr_specificity[rcc_tcr_specificity["Specificity"] == "NSUN5 p.Q18K"]
print(f"✅ Filtered for NSUN5 specificity: {rcc_tcr_specificity_filtered.shape[0]} rows retained.")

# Step 4: Create CDR3 to avidity mappings
cdr3a_to_avidity = dict(zip(rcc_tcr_specificity_filtered["CDR3A"], rcc_tcr_specificity_filtered["NSUN5_Avidity"]))
cdr3b_to_avidity = dict(zip(rcc_tcr_specificity_filtered["CDR3B"], rcc_tcr_specificity_filtered["NSUN5_Avidity"]))

# Step 5: Map to AnnData object
adata_merged_filtered_seurat.obs["NSUN5_avidity_TRA_cdr3"] = adata_merged_filtered_seurat.obs["TRA_cdr3"].map(cdr3a_to_avidity)
adata_merged_filtered_seurat.obs["NSUN5_avidity_TRB_cdr3"] = adata_merged_filtered_seurat.obs["TRB_cdr3"].map(cdr3b_to_avidity)

# Step 6: Summary statistics
print("\n🔍 Summary of Mapped Avidity Scores:")
print(f"✔ Mapped TRA CDR3s: {adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].notna().sum()}")
print(f"✔ Mapped TRB CDR3s: {adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].notna().sum()}")

print("\n🔍 Sample TRA mappings:")
print(adata_merged_filtered_seurat.obs[['TRA_cdr3', 'NSUN5_avidity_TRA_cdr3']].dropna().sample(5))

print("\n🔍 Sample TRB mappings:")
print(adata_merged_filtered_seurat.obs[['TRB_cdr3', 'NSUN5_avidity_TRB_cdr3']].dropna().sample(5))

# Step 7: Cell counts
num_cells = adata_merged_filtered_seurat.shape[0]
num_cells_rcc_108_2 = (adata_merged_filtered_seurat.obs['sample_id'] == "rcc_108_2").sum()
num_t_cells = (adata_merged_filtered_seurat.obs['manual_celltype_annotation'] == "T-cell").sum()

# Step 8: NSUN5-specific TCR counts
num_ns5_specific_tcr = ((adata_merged_filtered_seurat.obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR")).sum()
num_ns5_specific_tcr_avidity = (
    (adata_merged_filtered_seurat.obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR") &
    (
        adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].notna() |
        adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].notna()
    )
).sum()

# Step 9: Define strong vs weak
strong_mask = (
    (adata_merged_filtered_seurat.obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR") &
    (
        (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'] < 80)) |
        (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'] < 80))
    )
)

weak_mask = (
    (adata_merged_filtered_seurat.obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR") &
    (
        ((adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'] >= 80)) &
         (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].notna() & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'] >= 80))) |
        ((adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'].isna()) & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'] >= 80)) |
        ((adata_merged_filtered_seurat.obs['NSUN5_avidity_TRB_cdr3'].isna()) & (adata_merged_filtered_seurat.obs['NSUN5_avidity_TRA_cdr3'] >= 80))
    )
)

# Step 10: Spatial breakdown
num_ns5_avidity_strong_inside = (strong_mask & (adata_merged_filtered_seurat.obs['inside_region'] == True)).sum()
num_ns5_avidity_strong_outside = (strong_mask & (adata_merged_filtered_seurat.obs['inside_region'] == False)).sum()
num_ns5_avidity_weak_inside = (weak_mask & (adata_merged_filtered_seurat.obs['inside_region'] == True)).sum()
num_ns5_avidity_weak_outside = (weak_mask & (adata_merged_filtered_seurat.obs['inside_region'] == False)).sum()

# Step 11: Unclassified
unclassified_mask = (
    (adata_merged_filtered_seurat.obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR") &
    ~strong_mask & ~weak_mask
)
unclassified_count = unclassified_mask.sum()

# Step 12: Fisher’s Exact Test
fisher_table = [[num_ns5_avidity_strong_inside, num_ns5_avidity_weak_inside],
                [num_ns5_avidity_strong_outside, num_ns5_avidity_weak_outside]]
odds_ratio, p_value = fisher_exact(fisher_table)

# Step 13: Print results
print(f"\n📊 Total cells: {num_cells}")
print(f"📊 RCC 108_2 cells: {num_cells_rcc_108_2}")
print(f"📊 Total T cells: {num_t_cells}")
print(f"📊 NSUN5-specific TCRs: {num_ns5_specific_tcr}")
print(f"📊 NSUN5-specific TCRs with avidity: {num_ns5_specific_tcr_avidity}")
print(f"📊 Strong TCRs inside region: {num_ns5_avidity_strong_inside}")
print(f"📊 Strong TCRs outside region: {num_ns5_avidity_strong_outside}")
print(f"📊 Weak TCRs inside region: {num_ns5_avidity_weak_inside}")
print(f"📊 Weak TCRs outside region: {num_ns5_avidity_weak_outside}")
print(f"📊 Unclassified TCRs: {unclassified_count}")
print(f"\n🧪 Fisher’s Exact Test:")
print(f"✔ Odds Ratio: {odds_ratio:.3f}")
print(f"✔ P-value: {p_value:.3e}")

# Step 14: Print unclassified details
print("\n🔍 Unclassified NSUN5-specific TCRs:")
print(adata_merged_filtered_seurat.obs.loc[unclassified_mask, [
    'manual_celltype_annotation', 'NSUN5_MUT_TCR', 'NSUN5_avidity_TRA_cdr3', 'NSUN5_avidity_TRB_cdr3'
]])
