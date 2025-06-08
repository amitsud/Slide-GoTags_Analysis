# Author: Amit Sud
# Date: 1st May 2025
# Description: Compares SIINFEKL_WPRE mutation status derived from gene expression versus Nanopore sequencing.
#              Filters out NA entries, generates a 2x2 contingency table, and calculates the Phi coefficient
#              to measure the strength of association.
# Input:
#   - adata_filtered: AnnData object with the following columns in `.obs`:
#       • 'SIINFEKL_WPRE_status' (expression-based classification)
#       • 'SIINFEKL_WPRE_nanopore_status' (Nanopore sequencing classification)
# Output:
#   - A printed contingency table
#   - The Phi coefficient measuring concordance between the two status variables

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

# Extract relevant columns from the AnnData object
data = adata_filtered.obs
relevant_data = data[['SIINFEKL_WPRE_status', 'SIINFEKL_WPRE_nanopore_status']]

# Drop rows with missing (NA) values in either column (both pandas NA and string "NA")
filtered_data = relevant_data.dropna(subset=['SIINFEKL_WPRE_status', 'SIINFEKL_WPRE_nanopore_status'])
filtered_data = filtered_data[
    (filtered_data['SIINFEKL_WPRE_status'] != 'NA') &
    (filtered_data['SIINFEKL_WPRE_nanopore_status'] != 'NA')
]

# Create a 2x2 contingency table
contingency_table = pd.crosstab(
    filtered_data['SIINFEKL_WPRE_status'],
    filtered_data['SIINFEKL_WPRE_nanopore_status']
)
print("Contingency Table:")
print(contingency_table)

# Function to calculate the Phi coefficient from a 2x2 table
def calculate_phi(contingency_table):
    if contingency_table.shape == (2, 2):
        chi2, _, _, _ = chi2_contingency(contingency_table)
        n = contingency_table.to_numpy().sum()
        phi = np.sqrt(chi2 / n)
        return phi
    else:
        raise ValueError("Contingency table is not 2x2; Phi coefficient requires a 2x2 matrix.")

# Calculate and print the Phi coefficient
phi_value = calculate_phi(contingency_table)
print(f'The Phi coefficient is: {phi_value:.4f}')
