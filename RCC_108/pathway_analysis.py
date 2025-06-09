# Author: Amit Sud
# Date: 1st May 2025
# Description: This script extracts gene lists for the top enriched pathways using g:Profiler.
#              For each category in a pathway enrichment result dictionary, it retrieves genes
#              associated with the top N pathways by adjusted p-value.
# Input:
#   - Dictionary of pathway enrichment results (e.g., from a GSEA or g:Profiler run), where each value is a DataFrame
#     with columns including 'native', 'name', and 'FDR-adjusted'.
# Output:
#   - Dictionary mapping each selected pathway name to its associated gene set

from gprofiler import GProfiler
import pandas as pd

def extract_genes_from_pathway_ids(pathway_results, top_n=50, organism="hsapiens"):
    """
    Extract genes for the top N enriched pathways using pathway IDs.

    Parameters:
        pathway_results (dict): Output from pathway enrichment analysis.
        top_n (int): Number of top pathways to extract genes from.
        organism (str): Organism (default: "hsapiens" for human).

    Returns:
        dict: {pathway_name: set(genes)} mapping top pathways to their genes.
    """
    gp = GProfiler(return_dataframe=True)
    pathway_genes = {}

    for category, enrichment_df in pathway_results.items():
        top_pathways = enrichment_df.sort_values("FDR-adjusted").head(top_n)

        for _, row in top_pathways.iterrows():
            pathway_id = row['native']
            pathway_name = row['name']

            # Query g:Profiler to get associated genes
            pathway_info = gp.convert(organism=organism, query=[pathway_id])

            if not pathway_info.empty and 'name' in pathway_info.columns:
                genes = set(pathway_info['name'].dropna().tolist())
                pathway_genes[pathway_name] = genes
                print(f"✅ Extracted {len(genes)} genes for pathway: {pathway_name} ({pathway_id})")
            else:
                print(f"⚠️ No genes found for pathway: {pathway_name} ({pathway_id})")

    return pathway_genes


# Example usage (requires `pathway_results` to be defined in your session)
# top_pathway_genes = extract_genes_from_pathway_ids(pathway_results, top_n=50)
# for pathway, genes in list(top_pathway_genes.items())[:5]:
#     print(f"🔬 {pathway}: {len(genes)} genes ➡️ {list(genes)[:5]} ...")
