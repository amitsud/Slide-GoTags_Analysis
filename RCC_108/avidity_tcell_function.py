# Author: Amit Sud
# Date: 1st May 2025
# Description: This script calculates module scores (e.g., Exhaustion, Cytotoxicity, Chemokine)
#              in an AnnData object using provided gene sets, and compares them between groups
#              defined by NSUN5-specific TCR avidity and spatial localization. It generates violin
#              plots annotated with p-values and Cliff's delta effect sizes.
# Input:
#   - AnnData object (e.g., adata_merged_filtered_seurat) with NSUN5-specific metadata
#   - Tab-delimited file with gene modules (one column per module)
# Output:
#   - PDF violin plots comparing module scores between strong vs weak avidity groups,
#     and inside vs outside spatial regions

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from itertools import product
import os

def cliffs_delta(x, y):
    x, y = list(x), list(y)
    n_x, n_y = len(x), len(y)
    gt = sum(1 for xi, yi in product(x, y) if xi > yi)
    lt = sum(1 for xi, yi in product(x, y) if xi < yi)
    return (gt - lt) / (n_x * n_y)

def calculate_module_scores(adata, module_file_path, modules=["Exhaustion", "Cytotoxicity", "Chemokine/Chemokine_receptor"]):
    module_scores = pd.read_csv(module_file_path, delimiter='\t').dropna(how='all')
    for module in modules:
        if module not in module_scores.columns:
            raise ValueError(f"'{module}' column not found.")
        gene_list = module_scores[module].dropna().tolist()
        valid_genes = [gene for gene in gene_list if gene in adata.var_names]
        if not valid_genes:
            raise ValueError(f"No valid genes in AnnData for module '{module}'.")
        sc.tl.score_genes(adata, gene_list=valid_genes, score_name=module)
        print(f"✅ Calculated module score for '{module}' using {len(valid_genes)} genes.")

def compare_and_plot_separately(adata, module_file, save_path="module_score_plots"):
    os.makedirs(save_path, exist_ok=True)
    calculate_module_scores(adata, module_file)

    obs = adata.obs.copy()

    # Define Avidity
    strong = (
        (obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR") &
        (
            (obs['NSUN5_avidity_TRA_cdr3'].notna() & (obs['NSUN5_avidity_TRA_cdr3'] < 80)) |
            (obs['NSUN5_avidity_TRB_cdr3'].notna() & (obs['NSUN5_avidity_TRB_cdr3'] < 80))
        )
    )
    weak = (
        (obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR") &
        (
            ((obs['NSUN5_avidity_TRA_cdr3'].notna() & (obs['NSUN5_avidity_TRA_cdr3'] >= 80)) &
             (obs['NSUN5_avidity_TRB_cdr3'].notna() & (obs['NSUN5_avidity_TRB_cdr3'] >= 80))) |
            ((obs['NSUN5_avidity_TRA_cdr3'].isna()) & (obs['NSUN5_avidity_TRB_cdr3'] >= 80)) |
            ((obs['NSUN5_avidity_TRB_cdr3'].isna()) & (obs['NSUN5_avidity_TRA_cdr3'] >= 80))
        )
    )
    obs['Avidity'] = np.nan
    obs.loc[strong, 'Avidity'] = "Strong avidity"
    obs.loc[weak, 'Avidity'] = "Weak avidity"

    obs = obs[obs['NSUN5_MUT_TCR'] == "NSUN5_specific_TCR"].copy()
    obs = obs.dropna(subset=["Avidity", "inside_region", "Exhaustion", "Cytotoxicity", "Chemokine/Chemokine_receptor"])

    comparisons = {
        "High_vs_Low_Avidity": obs.copy(),
        "Inside_vs_Outside_Region": obs.copy(),
        "High_vs_Low_Avidity_Inside": obs[obs['inside_region'] == True].copy(),
        "High_vs_Low_Avidity_Outside": obs[obs['inside_region'] == False].copy()
    }

    # 🎨 Color palette
    palette = {
        "Strong avidity": "#ec3037",
        "Weak avidity": "#c14eac",
        "Inside": "#ec3037",
        "Outside": "#3f5ef8"
    }

    for module in ["Exhaustion", "Cytotoxicity", "Chemokine/Chemokine_receptor"]:
        for label, data in comparisons.items():
            plt.figure(figsize=(5, 5))

            if label == "Inside_vs_Outside_Region":
                data['Region'] = np.where(data['inside_region'], "Inside", "Outside")
                if data['Region'].nunique() < 2:
                    print(f"⚠️ Skipping {label} for {module} (not enough region groups)")
                    continue
                x = "Region"
                groups = ["Inside", "Outside"]
                colors = [palette["Inside"], palette["Outside"]]
            else:
                if data['Avidity'].nunique() < 2:
                    print(f"⚠️ Skipping {label} for {module} (not enough avidity groups)")
                    continue
                x = "Avidity"
                groups = ["Strong avidity", "Weak avidity"]
                colors = [palette["Strong avidity"], palette["Weak avidity"]]

            group1 = data[data[x] == groups[0]][module]
            group2 = data[data[x] == groups[1]][module]

            if group1.empty or group2.empty:
                print(f"⚠️ Skipping {label} for {module} (empty group)")
                continue

            ax = sns.violinplot(
                data=data,
                x=x,
                y=module,
                palette=colors,
                order=groups,
                inner=None,
                linewidth=2,
                cut=0,
                width=0.7
            )

            sns.stripplot(
                data=data,
                x=x,
                y=module,
                color="black",
                size=5,
                jitter=False,
                alpha=0.6,
                order=groups
            )

            ax.grid(False, axis='y')

            # Stats
            stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
            delta = cliffs_delta(group1, group2)

            plt.title(f"{module} — {label.replace('_', ' ')}", fontsize=11)
            plt.ylabel(f"{module} Module Score")
            plt.xlabel("")
            plt.xticks(fontsize=10)

            text = f"p = {p:.3f}\nCliff’s δ = {delta:.2f}"
            plt.text(0.5, max(data[module]), text, ha='center', va='top', fontsize=9)

            sns.despine()
            plt.tight_layout()

            fname = f"{save_path}/{module.replace('/', '_')}_{label}.pdf"
            plt.savefig(fname, format='pdf', dpi=1200)
            print(f"📄 Saved: {fname}")
            plt.show()

# Run
module_file = "path/to/module_scores.txt"
compare_and_plot_separately(adata_merged_filtered_seurat, module_file)
