"""UMAP visualization with cluster annotations."""
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

# Snakemake parameters
input_h5ad = str(snakemake.input.h5ad)
output_umap = str(snakemake.output.umap)
cluster_key = snakemake.params.cluster_key
log_file = str(snakemake.log[0])

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Loading {input_h5ad}...")
adata = sc.read_h5ad(input_h5ad)

# Set figure parameters
sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(8, 8))

# Plot UMAP
print(f"Plotting UMAP colored by {cluster_key}...")
fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(adata, color=cluster_key, legend_loc='on data', 
           title=f'UMAP - {cluster_key} clusters', ax=ax, show=False)
plt.tight_layout()
plt.savefig(output_umap, dpi=150, bbox_inches='tight')
plt.close()

print(f"Saved UMAP to {output_umap}")