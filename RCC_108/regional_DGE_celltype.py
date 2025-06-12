# Author: Amit Sud
# Description: This script performs differential gene expression (DEG) analysis across spatial regions
#              within each cell type in a single-cell dataset. It highlights genes from selected pathways
#              and visualizes results as volcano plots with annotations.
# Input:
#   - AnnData object (`adata_merged_filtered_seurat`) with:
#       - 'sample_id', 'manual_celltype_annotation', and 'inside_region' in `.obs`
#   - Dictionary of pathway genes (e.g., from EnrichR or MSigDB)
#   - List of selected pathways (e.g., ["Antigen processing and presentation", "Interferon gamma signaling"])
# Output:
#   - Volcano plots (PDFs) per cell type showing DEG results with colored highlights for selected pathways

from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os

def perform_deg_and_plot(
    adata, sample_id, group_column, region_column, 
    fdr_threshold=0.05, min_cells_expressing=5, 
    min_mean_expression=0.0, logfc_clip=None, pseudocount=1e-6,
    pathway_genes=None, selected_pathways=None, output_dir="./"
):
    os.makedirs(output_dir, exist_ok=True)

    if pathway_genes and selected_pathways:
        pathway_genes = {k: v for k, v in pathway_genes.items() if k in selected_pathways}
        print(f"🔎 Using pathways: {list(pathway_genes.keys())}")
    else:
        pathway_genes = {}

    pathway_sets = [pathway_genes.get(p, set()) for p in selected_pathways] if selected_pathways else []
    overlap_genes = set.intersection(*pathway_sets) if len(pathway_sets) > 1 else set()

    sample_data = adata[adata.obs['sample_id'] == sample_id].copy()
    categories = sample_data.obs[group_column].unique()
    deg_results = {}

    for category in categories:
        category_data = sample_data[sample_data.obs[group_column] == category]
        inside = category_data[category_data.obs[region_column] == True].to_df()
        outside = category_data[category_data.obs[region_column] == False].to_df()

        print(f"\n🔹 Category: {category}")
        print(f"📌 Cells inside: {inside.shape[0]}, outside: {outside.shape[0]}")

        if inside.shape[0] < 2 or outside.shape[0] < 2:
            print(f"⚠ Skipping {category}: insufficient cells.")
            continue

        genes = sample_data.var.index
        log_fold_changes, p_values, retained_genes = [], [], []

        if isinstance(min_cells_expressing, float):
            inside_threshold = int(np.ceil(min_cells_expressing * inside.shape[0]))
            outside_threshold = int(np.ceil(min_cells_expressing * outside.shape[0]))
        elif isinstance(min_cells_expressing, int):
            inside_threshold = outside_threshold = min_cells_expressing
        else:
            raise ValueError("min_cells_expressing must be integer ≥1 or float between 0 and 1.")

        for gene in genes:
            inside_expr = inside[gene]
            outside_expr = outside[gene]
            inside_mean = inside_expr.mean()
            outside_mean = outside_expr.mean()
            inside_detected = (inside_expr > 0).sum()
            outside_detected = (outside_expr > 0).sum()

            if (inside_detected < inside_threshold and outside_detected < outside_threshold):
                continue
            if (inside_mean < min_mean_expression and outside_mean < min_mean_expression):
                continue

            log_fc = np.log2((inside_mean + pseudocount) / (outside_mean + pseudocount))
            if logfc_clip:
                log_fc = np.clip(log_fc, -logfc_clip, logfc_clip)

            try:
                _, p_val = ttest_ind(inside_expr, outside_expr, equal_var=False, nan_policy='omit')
            except ValueError:
                p_val = np.nan

            retained_genes.append(gene)
            log_fold_changes.append(log_fc)
            p_values.append(p_val)

        if not retained_genes:
            print("⚠ No genes passed filtering.")
            continue

        fdrs = np.full(len(p_values), np.nan)
        valid_p = np.array(p_values)[~np.isnan(p_values)]
        if valid_p.size > 0:
            fdrs[~np.isnan(p_values)] = multipletests(valid_p, method='fdr_bh')[1]

        results_df = pd.DataFrame({
            "Gene": retained_genes,
            "Log2FC": log_fold_changes,
            "P-Value": p_values,
            "FDR": fdrs
        }).sort_values("FDR", na_position='last')

        deg_results[category] = results_df
        sig_genes = results_df[results_df["FDR"] < fdr_threshold]
        print(f"🔍 Significant genes (FDR < {fdr_threshold}): {len(sig_genes)}")

        fig, ax = plt.subplots(figsize=(14, 8))
        ax.scatter(results_df['Log2FC'], -np.log10(results_df['P-Value']),
                   c='grey', alpha=0.5, label='Non-significant', s=20)

        pathway_colors = {
            selected_pathways[0]: 'blue',
            selected_pathways[1]: 'green',
            "Both": 'teal',
            "Significant": 'red'
        }

        for label, color in pathway_colors.items():
            if label == "Both":
                subset = results_df[results_df["Gene"].isin(overlap_genes)]
            elif label == "Significant":
                subset = sig_genes[~sig_genes["Gene"].isin(overlap_genes)]
            else:
                subset = results_df[
                    results_df["Gene"].isin(pathway_genes.get(label, set())) &
                    ~results_df["Gene"].isin(overlap_genes)
                ]

            ax.scatter(subset["Log2FC"], -np.log10(subset["P-Value"]),
                       c=color, label=label, s=40, edgecolors='black', alpha=0.8)

        texts = [
            ax.text(row['Log2FC'], -np.log10(row['P-Value']) + 0.2, row['Gene'],
                    fontsize=8, ha='center', va='bottom')
            for _, row in sig_genes.iterrows()
        ]
        adjust_text(texts, arrowprops=dict(arrowstyle='-', lw=0.5, color='black'))

        ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1.2, label='p = 0.05')
        max_logfc = max(abs(results_df['Log2FC'].min()), results_df['Log2FC'].max())
        ax.set_xlim(-max_logfc - 0.5, max_logfc + 0.5)
        ax.set_title(f"Volcano Plot: {category} ({sample_id})", fontsize=14)
        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-Log10 P-Value")
        ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
        plt.tight_layout(rect=[0, 0, 0.85, 1])

        pdf_path = os.path.join(output_dir, f"volcano_plot_pathways_{category}_{sample_id}.pdf")
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"✅ Saved: {pdf_path}")

    return deg_results


# Example usage
selected_pathways = ["Antigen processing and presentation", "Interferon gamma signaling"]

deg_results = perform_deg_and_plot(
    adata=adata_merged_filtered_seurat,
    sample_id='rcc_108_2',
    group_column='manual_celltype_annotation',
    region_column='inside_region',
    fdr_threshold=0.05,
    min_cells_expressing=0.2,
    min_mean_expression=0.1,
    pathway_genes=top_pathway_genes,
    selected_pathways=selected_pathways,
    output_dir="./"
)
