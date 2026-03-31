"""Normalization and HVG selection using Scanpy."""
import scanpy as sc
import sys

# Snakemake parameters - extract actual string paths
input_h5ad = str(snakemake.input.h5ad)
output_h5ad = str(snakemake.output[0])
target_sum = snakemake.params.target_sum
n_top_genes = snakemake.params.n_top_genes
log_file = str(snakemake.log[0])

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Loading {input_h5ad}...")
adata = sc.read_h5ad(input_h5ad)
print(f"Shape: {adata.shape}")

# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize
print(f"Normalizing to target_sum={target_sum}...")
sc.pp.normalize_total(adata, target_sum=target_sum)

# Log transform
print("Log1p transformation...")
sc.pp.log1p(adata)

# Highly variable genes
print(f"Selecting {n_top_genes} highly variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=False)
print(f"HVGs selected: {adata.var.highly_variable.sum()}")

# Save
adata.write_h5ad(output_h5ad)
print(f"Saved to {output_h5ad}")