# Author: Amit Sud
# Description: Visualizes the frequency of TRA and TRB clonotypes in T cells.
#              Highlights specific clonotypes of interest, and presents counts in side-by-side bar charts.
# Input:
#   - adata_filtered: AnnData object with `.obs` columns:
#       • 'TRA_filtered' and 'TRB_filtered' (clonotype group identifiers)
# Output:
#   - Side-by-side bar charts of TRA and TRB clonotype frequencies
#   - (Optional) A PDF file "TRA_TRB_Clonotypes_Split_Grid_Highlighted.pdf"

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Extract clonotype groupings for TRA and TRB
TRA_filtered = adata_filtered.obs['TRA_filtered']
TRB_filtered = adata_filtered.obs['TRB_filtered']

# Combine TRA and TRB clonotype frequencies into a single dataframe
cdr3_combined = pd.DataFrame({
    'TRA': TRA_filtered.value_counts(),
    'TRB': TRB_filtered.value_counts()
}).fillna(0).astype(int)

# Add total count column for sorting
cdr3_combined['Total'] = cdr3_combined['TRA'] + cdr3_combined['TRB']

# Sort by TRA then TRB frequency (descending)
cdr3_combined = cdr3_combined.sort_values(by=['TRA', 'TRB'], ascending=[False, False])

# Remove any placeholder or uninformative entries
cdr3_combined = cdr3_combined.loc[~cdr3_combined.index.isin(['Not genotyped'])]

# Define specific clonotypes to highlight
special_values = ['CAASDNYQLIW', 'CASSRANYEQYF']

# Separate TRA and TRB for individual plotting
cdr3_tra = cdr3_combined[['TRA']].copy()
cdr3_tra = cdr3_tra[cdr3_tra['TRA'] > 0]

cdr3_trb = cdr3_combined[['TRB']].copy()
cdr3_trb = cdr3_trb[cdr3_trb['TRB'] > 0]

# Set colors: red for special clonotypes, light red for others
bar_colors_tra = ['#FF0000' if clonotype in special_values else '#FFCCCC' for clonotype in cdr3_tra.index]
bar_colors_trb = ['#FF0000' if clonotype in special_values else '#FFCCCC' for clonotype in cdr3_trb.index]

# Adjust figure size based on number of bars
bar_width = 0.8
fig_width_tra = len(cdr3_tra) * bar_width * 0.6
fig_width_trb = len(cdr3_trb) * bar_width * 0.6
fig = plt.figure(figsize=(fig_width_tra + fig_width_trb, 8))
gs = fig.add_gridspec(1, 2, width_ratios=[fig_width_tra, fig_width_trb], wspace=0.3)

# --- Plot TRA clonotypes ---
ax1 = fig.add_subplot(gs[0])
ax1.bar(range(len(cdr3_tra)), cdr3_tra['TRA'], color=bar_colors_tra, edgecolor='black', width=bar_width)
ax1.set_title('TRA Clonotypes', fontsize=18)
ax1.set_xlabel('TRA Clonotype', fontsize=14)
ax1.set_ylabel('Number of T cells', fontsize=14)
ax1.set_xticks(range(len(cdr3_tra.index)))
ax1.set_xticklabels(cdr3_tra.index, rotation=90, fontsize=10)
ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
ax1.grid(False)

# --- Plot TRB clonotypes ---
ax2 = fig.add_subplot(gs[1])
ax2.bar(range(len(cdr3_trb)), cdr3_trb['TRB'], color=bar_colors_trb, edgecolor='black', width=bar_width)
ax2.set_title('TRB Clonotypes', fontsize=18)
ax2.set_xlabel('TRB Clonotype', fontsize=14)
ax2.set_ylabel('Number of T cells', fontsize=14)
ax2.set_xticks(range(len(cdr3_trb.index)))
ax2.set_xticklabels(cdr3_trb.index, rotation=90, fontsize=10)
ax2.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
ax2.grid(False)

# --- Synchronize y-axis limits ---
y_max = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
ax1.set_ylim(0, y_max)
ax2.set_ylim(0, y_max)

# --- Add legend to explain colors ---
legend_patches = [
    Patch(facecolor='#FF0000', edgecolor='black', label='Special Clonotype', linewidth=1),
    Patch(facecolor='#FFCCCC', edgecolor='black', label='Other Clonotype', linewidth=1)
]
fig.legend(
    handles=legend_patches,
    title='Clonotype Highlights',
    fontsize=12,
    title_fontsize=12,
    loc='upper center',
    ncol=2,
    bbox_to_anchor=(0.5, 1.05),
    frameon=True
)

# Adjust layout to fit plots and legend
plt.tight_layout(rect=[0, 0, 1, 0.95])

# --- Save and show plot ---
# Uncomment below to save as PDF
# plt.savefig("TRA_TRB_Clonotypes_Split_Grid_Highlighted.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.show()
