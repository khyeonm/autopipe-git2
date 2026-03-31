"""Marker genes dotplot visualization."""
import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

# Snakemake parameters
input_h5ad = str(snakemake.input.h5ad)
input_markers = str(snakemake.input.markers)
output_dotplot = str(snakemake.output.dotplot)
n_genes = snakemake.params.n_genes
cluster_key = snakemake.params.cluster_key
log_file = str(snakemake.log[0])

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Loading {input_h5ad}...")
adata = sc.read_h5ad(input_h5ad)

print(f"Loading markers from {input_markers}...")
markers_df = pd.read_csv(input_markers)

# Get top n genes per cluster (unique, ordered by score)
top_genes = markers_df[markers_df['rank'] <= n_genes]['gene'].unique().tolist()
print(f"Plotting dotplot for {len(top_genes)} genes...")

# Set figure parameters
sc.settings.set_figure_params(dpi=150, frameon=False)

# Plot dotplot
sc.pl.rank_genes_groups_dotplot(
    adata,
    n_genes=n_genes,
    groupby=cluster_key,
    standard_scale='var',
    show=False
)
plt.savefig(output_dotplot, dpi=150, bbox_inches='tight')
plt.close()

print(f"Saved dotplot to {output_dotplot}")