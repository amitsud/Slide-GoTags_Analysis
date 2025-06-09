# Author: Amit Sud  
# Date: 1st May 2025  
# Description: This script generates spatial plots for each sample, 
#              coloring cells by their NSUN5_MUT_TCR classification.
# Input:
#   - AnnData object: adata_merged_filtered_seurat
#     Required columns in `.obs`: 'sample_id', 'x', 'y', 'NSUN5_MUT_TCR'
# Output:
#   - One PDF spatial plot per sample showing cell classifications

import matplotlib.pyplot as plt

# --- Define Color Mapping and Plot Order ---
color_mapping = {
    'NSUN5_specific_TCR': '#FF0000',   # Red
    'Tumor_specific_TCR': '#F8C47D',   # Light Orange
    'Other_clonotype_TCR': '#FFCCCC',  # Light Red
    'NSUN5_MUT': '#333399',            # Dark Blue
    'NSUN5_WT': '#A5CEE3',             # Light Blue
    'NA': '#D9D9D9'                    # Grey for NA
}

plot_order = ['NA', 'NSUN5_WT', 'Other_clonotype_TCR', 'Tumor_specific_TCR', 'NSUN5_specific_TCR', 'NSUN5_MUT']

# --- Generate Plots by Sample ---
if 'sample_id' in adata_merged_filtered_seurat.obs:
    unique_samples = adata_merged_filtered_seurat.obs['sample_id'].unique()

    for sample_id in unique_samples:
        sample_data = adata_merged_filtered_seurat[adata_merged_filtered_seurat.obs['sample_id'] == sample_id]

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_aspect('equal', adjustable='datalim')

        # Plot NA background
        na_subset = sample_data[sample_data.obs['NSUN5_MUT_TCR'].isna()]
        ax.scatter(
            na_subset.obs['x'],
            na_subset.obs['y'],
            c=color_mapping['NA'],
            edgecolor='none',
            s=50,
            alpha=0.9,
            label=f"NA (n={na_subset.shape[0]})"
        )

        # Plot defined categories
        for category in plot_order:
            if category != 'NA' and category in color_mapping:
                subset = sample_data[sample_data.obs['NSUN5_MUT_TCR'] == category]
                ax.scatter(
                    subset.obs['x'],
                    subset.obs['y'],
                    c=color_mapping[category],
                    edgecolor='black',
                    s=50,
                    alpha=0.9,
                    label=f"{category} (n={subset.shape[0]})",
                    linewidths=0.5
                )

        # Add scale bars
        micron_to_pixel_ratio = 1  # Adjust if needed
        size_microns = 500
        size_pixels = size_microns * micron_to_pixel_ratio

        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        x0 = x_min + (x_max - x_min) * 0.02
        y0 = y_min + (y_max - y_min) * 0.02

        ax.plot([x0, x0 + size_pixels], [y0, y0], 'k-', linewidth=2)
        ax.text(x0 + size_pixels / 2, y0 - (y_max - y_min) * 0.015,
                f'{size_microns} µm', ha='center', va='top', fontsize=10)
        ax.plot([x0, x0], [y0, y0 + size_pixels], 'k-', linewidth=2)
        ax.text(x0 - (x_max - x_min) * 0.01, y0 + size_pixels / 2,
                f'{size_microns} µm', ha='right', va='center', fontsize=10, rotation=90)

        # Final plot formatting
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, title='NSUN5_MUT_TCR')
        ax.set_title(f'Spatial Plot for Sample {sample_id}', fontsize=16, pad=20)
        ax.set_xlabel('Spatial X', fontsize=12, labelpad=10)
        ax.set_ylabel('Spatial Y', fontsize=12, labelpad=10)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)

        for spine in ax.spines.values():
            spine.set_linewidth(2)

        plt.tight_layout()

        # Save plot to generic path
        output_file = f"path/to/output/Spatial_Plot_NSUN5_MUT_TCR_Sample_{sample_id}.pdf"
        # plt.savefig(output_file, dpi=1200, bbox_inches='tight', format='pdf')
        # plt.show()
        plt.close(fig)
        print(f"Spatial plot saved as {output_file}")

else:
    print("Error: Required column 'sample_id' not found in adata_merged_filtered_seurat.obs.")
