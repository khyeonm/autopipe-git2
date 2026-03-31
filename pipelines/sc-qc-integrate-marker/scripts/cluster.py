"""PCA, UMAP, and Leiden clustering using Scanpy."""
import scanpy as sc
import sys

# Snakemake parameters
input_h5ad = str(snakemake.input.h5ad)
output_h5ad = str(snakemake.output[0])
n_pcs = snakemake.params.n_pcs
n_neighbors = snakemake.params.n_neighbors
resolution = snakemake.params.resolution
log_file = str(snakemake.log[0])

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Loading {input_h5ad}...")
adata = sc.read_h5ad(input_h5ad)
print(f"Shape: {adata.shape}")

# Scale data (on HVGs)
print("Scaling data...")
sc.pp.scale(adata, max_value=10)

# PCA
print(f"Running PCA with {n_pcs} components...")
sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
print(f"Variance explained (top 10): {adata.uns['pca']['variance_ratio'][:10].sum():.2%}")

# Neighbors graph
print(f"Computing neighbors with n_neighbors={n_neighbors}...")
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

# UMAP
print("Computing UMAP...")
sc.tl.umap(adata)

# Leiden clustering
print(f"Running Leiden clustering with resolution={resolution}...")
sc.tl.leiden(adata, resolution=resolution)
print(f"Found {adata.obs['leiden'].nunique()} clusters")

# Save
adata.write_h5ad(output_h5ad)
print(f"Saved to {output_h5ad}")