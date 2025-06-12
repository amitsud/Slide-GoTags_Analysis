# Author: Amit Sud
# Description: This script extracts TCR clonotype information from spatial single-cell data,
#              computes relative frequencies in two regions (defined by a polygon mask), and
#              visualizes these as paired stacked bar plots with shaded connections.
# Input:
#   - AnnData object (adata_merged_filtered_seurat) with:
#       - 'sample_id', 'manual_celltype_annotation', 'inside_region',
#         'TRA_cdr3', 'TRB_cdr3', 'TRA_filtered', 'TRB_filtered' in obs
# Output:
#   - PDF files with stacked bar plots of clonotype distributions for TRA and TRB chains

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

# ✅ Color mapping based on filtered assignments
color_map = {
    'Other clonotype': '#FFCCCC',
    'Tumor_specific_TRA': '#F8C47D',
    'Tumor_specific_TRB': '#F8C47D',
    'NSUN5_specific_TRA': '#FF0000',
    'NSUN5_specific_TRB': '#FF0000'
}

# ✅ Generate lighter color for shading
def lighter_color(hex_color, factor=0.8):
    rgb = mcolors.hex2color(hex_color)
    lighter_rgb = [1 - (1 - c) * factor for c in rgb]
    return mcolors.to_hex(lighter_rgb)

# ✅ Group order for stacking
group_order = ['NSUN5_specific_TRA', 'NSUN5_specific_TRB', 'Tumor_specific_TRA', 'Tumor_specific_TRB', 'Other clonotype']

# ✅ Extract TCR sequences using inside_region and filter for T-cells
def extract_tcr_within_polygon(adata, sample_id):
    sample_data = adata[adata.obs['sample_id'] == sample_id].copy()
    sample_data.obs['region_assignment'] = sample_data.obs['inside_region'].map({True: 'Region 1', False: 'Region 2'})
    t_cell_data = sample_data[sample_data.obs['manual_celltype_annotation'] == 'T-cell']

    tcr_df = t_cell_data.obs[['region_assignment', 'TRA_cdr3', 'TRB_cdr3', 'TRA_filtered', 'TRB_filtered']].copy()
    for col in ['TRA_cdr3', 'TRB_cdr3', 'TRA_filtered', 'TRB_filtered']:
        tcr_df[col] = tcr_df[col].astype(str).replace('nan', 'None')
    return tcr_df

# ✅ Calculate clonotype frequencies
def calculate_clonotype_frequencies(tcr_df, chain, filter_col):
    region1_df = tcr_df[(tcr_df['region_assignment'] == 'Region 1') & (tcr_df[chain] != 'None')]
    region1_counts = region1_df.groupby([chain, filter_col]).size().reset_index(name='count')
    region1_counts['region'] = 'Region 1'
    region1_counts['relative_frequency'] = region1_counts['count'] / region1_counts['count'].sum()

    target_clonotypes = region1_df[chain].unique()
    region2_df = tcr_df[(tcr_df['region_assignment'] == 'Region 2') & (tcr_df[chain].isin(target_clonotypes))]
    total_region2_clonotypes = tcr_df[(tcr_df['region_assignment'] == 'Region 2') & (tcr_df[chain] != 'None')].shape[0]
    region2_counts = region2_df.groupby([chain, filter_col]).size().reset_index(name='count')
    region2_counts['region'] = 'Region 2'
    region2_counts['relative_frequency'] = region2_counts['count'] / total_region2_clonotypes

    freq_df = pd.concat([region1_counts, region2_counts], ignore_index=True)
    freq_df[filter_col] = pd.Categorical(freq_df[filter_col], categories=group_order, ordered=True)
    freq_df = freq_df.sort_values(by=[filter_col, chain, 'relative_frequency'], ascending=[True, True, False])
    return freq_df

# ✅ Plot stacked bar graph
def plot_stacked_bar(freq_df, chain, filter_col, save_pdf=False, pdf_filename="stacked_bar_plot.pdf"):
    pivot_df = freq_df.pivot_table(index=[chain, filter_col], columns='region', values='relative_frequency', fill_value=0).reset_index()
    pivot_df[filter_col] = pd.Categorical(pivot_df[filter_col], categories=group_order, ordered=True)
    pivot_df = pivot_df.sort_values(by=[filter_col, 'Region 1'], ascending=[True, False])

    fig, ax = plt.subplots(figsize=(6, 8))
    bar_positions = {'Region 1': -0.03, 'Region 2': 0.07}
    bar_width = 0.08
    outer_line_width = 2.0
    inner_line_width = 0.7

    bottom_r1, bottom_r2 = 0, 0
    lines_data = []

    for _, row in pivot_df.iterrows():
        clonotype = row[chain]
        group = row[filter_col]
        base_color = color_map.get(group, '#D3D3D3')
        fill_color = lighter_color(base_color, factor=0.6)

        p1, p2 = row['Region 1'], row['Region 2']
        ax.bar(bar_positions['Region 1'], p1, bottom=bottom_r1, width=bar_width, color=base_color,
               edgecolor='black', linewidth=inner_line_width, zorder=3)
        ax.bar(bar_positions['Region 2'], p2, bottom=bottom_r2, width=bar_width, color=base_color,
               edgecolor='black', linewidth=inner_line_width, zorder=3)

        x1_right = bar_positions['Region 1'] + bar_width / 2
        x2_left = bar_positions['Region 2'] - bar_width / 2
        y1_top, y2_top = bottom_r1 + p1, bottom_r2 + p2

        polygon_verts = [
            (x1_right, bottom_r1), (x1_right, y1_top),
            (x2_left, y2_top), (x2_left, bottom_r2), (x1_right, bottom_r1)
        ]
        polygon = mpatches.Polygon(polygon_verts, closed=True, facecolor=fill_color, alpha=0.9, edgecolor='none', zorder=1)
        ax.add_patch(polygon)

        lines_data.append(((x1_right, y1_top), (x2_left, y2_top)))
        bottom_r1 += p1
        bottom_r2 += p2

    for region, bottom in zip(['Region 1', 'Region 2'], [bottom_r1, bottom_r2]):
        pos = bar_positions[region]
        ax.plot([pos - bar_width / 2, pos + bar_width / 2], [0, 0], color='black', linewidth=outer_line_width, zorder=4)
        ax.plot([pos - bar_width / 2, pos + bar_width / 2], [bottom, bottom], color='black', linewidth=outer_line_width, zorder=4)
        ax.plot([pos - bar_width / 2, pos - bar_width / 2], [0, bottom], color='black', linewidth=outer_line_width, zorder=4)
        ax.plot([pos + bar_width / 2, pos + bar_width / 2], [0, bottom], color='black', linewidth=outer_line_width, zorder=4)

    for (x1, y1), (x2, y2) in lines_data:
        ax.plot([x1, x2], [y1, y2], color='black', linewidth=0.8, alpha=0.9, zorder=5)

    ax.set_xticks([bar_positions['Region 1'], bar_positions['Region 2']])
    ax.set_xticklabels(['Region 1', 'Region 2'], fontsize=12)
    ax.set_ylabel("Relative Frequency", fontsize=12)
    ax.set_title(f"{chain} Clonotype Distribution", fontsize=14)
    ax.set_ylim(0, 1)
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(False)
    plt.tight_layout()

    if save_pdf:
        # plt.savefig(pdf_filename, dpi=1200, bbox_inches='tight')
        print(f"✅ Plot saved as: {pdf_filename}")

# ✅ Example usage
tcr_df = extract_tcr_within_polygon(
    adata=adata_merged_filtered_seurat,
    sample_id='sample_id_placeholder'
)

if tcr_df is not None:
    plot_configs = [
        ('TRA_cdr3', 'TRA_filtered', "stacked_bar_TRA_cdr3_clonotype_distribution.pdf"),
        ('TRB_cdr3', 'TRB_filtered', "stacked_bar_TRB_cdr3_clonotype_distribution.pdf")
    ]

    for chain, filter_col, filename in plot_configs:
        freq_df = calculate_clonotype_frequencies(tcr_df, chain, filter_col)
        plot_stacked_bar(freq_df, chain=chain, filter_col=filter_col, save_pdf=True, pdf_filename=filename)
        
