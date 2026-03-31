"""Automatic marker gene detection using Scanpy."""
import scanpy as sc
import pandas as pd
import sys

# Snakemake parameters
input_h5ad = str(snakemake.input.h5ad)
output_markers = str(snakemake.output.markers)
output_h5ad = str(snakemake.output.h5ad)
n_genes = snakemake.params.n_genes
method = snakemake.params.method
cluster_key = snakemake.params.cluster_key
log_file = str(snakemake.log[0])

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Loading {input_h5ad}...")
adata = sc.read_h5ad(input_h5ad)
print(f"Shape: {adata.shape}")
print(f"Clusters in '{cluster_key}': {adata.obs[cluster_key].nunique()}")

# Find marker genes
print(f"Finding marker genes using {method} test...")
sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method, n_genes=n_genes)

# Extract results to DataFrame
print("Extracting marker results...")
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names

markers_list = []
for group in groups:
    for i in range(n_genes):
        markers_list.append({
            'cluster': group,
            'rank': i + 1,
            'gene': result['names'][group][i],
            'score': result['scores'][group][i],
            'logfoldchange': result['logfoldchanges'][group][i],
            'pval': result['pvals'][group][i],
            'pval_adj': result['pvals_adj'][group][i]
        })

markers_df = pd.DataFrame(markers_list)
markers_df.to_csv(output_markers, index=False)
print(f"Saved markers to {output_markers}")

# Save annotated h5ad
adata.write_h5ad(output_h5ad)
print(f"Saved annotated data to {output_h5ad}")