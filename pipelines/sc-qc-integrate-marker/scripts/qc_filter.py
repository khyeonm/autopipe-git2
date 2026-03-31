"""QC filtering for single-cell RNA-seq data using Scanpy."""
import scanpy as sc
import sys

# Snakemake parameters - extract actual string paths
input_h5ad = str(snakemake.input.h5ad)
output_h5ad = str(snakemake.output[0])
min_genes = snakemake.params.min_genes
min_cells = snakemake.params.min_cells
max_genes = snakemake.params.max_genes
max_mt_pct = snakemake.params.max_mt_pct
log_file = str(snakemake.log[0])

# Redirect output to log
sys.stdout = open(log_file, 'w')
sys.stderr = sys.stdout

print(f"Loading {input_h5ad}...")
adata = sc.read_h5ad(input_h5ad)
print(f"Initial shape: {adata.shape}")

# Basic QC metrics
print("Calculating QC metrics...")
adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter cells
print(f"Filtering cells: min_genes={min_genes}, max_genes={max_genes}, max_mt_pct={max_mt_pct}")
sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
adata = adata[adata.obs.pct_counts_mt < max_mt_pct, :]

# Filter genes
print(f"Filtering genes: min_cells={min_cells}")
sc.pp.filter_genes(adata, min_cells=min_cells)

print(f"After QC filtering: {adata.shape}")
adata.write_h5ad(output_h5ad)
print(f"Saved to {output_h5ad}")