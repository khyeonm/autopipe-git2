# sc-qc-integrate-marker

Hybrid single-cell RNA-seq analysis pipeline: QC/normalization with Scanpy → Integration/clustering with nf-core/scdownstream → Automatic marker gene detection and visualization.

## Pipeline Steps

| Step | Tool | Description |
|------|------|-------------|
| 1 | Scanpy | QC filtering + Normalization + HVG selection |
| 2 | nf-core/scdownstream | Dimensionality reduction + Clustering |
| 3 | Scanpy | Marker gene detection + Visualization |

## Input

- **h5ad file**: Raw single-cell count matrix in AnnData format

## Outputs

| File | Description |
|------|-------------|
| `qc_normalized.h5ad` | QC-filtered and normalized data (Step 1) |
| `clustered.h5ad` | Clustered data with embeddings (Step 2) |
| `marker_analysis/cluster_markers.csv` | Marker genes per cluster with statistics |
| `marker_analysis/annotated.h5ad` | Final annotated data |
| `marker_analysis/umap_clusters.png` | UMAP colored by cluster |
| `marker_analysis/top_markers_heatmap.png` | Heatmap of top markers |
| `marker_analysis/top_markers_dotplot.png` | Dotplot of top markers |

## Configuration

Edit `config.yaml` to customize:

### QC Parameters
- `min_genes`: Minimum genes per cell (default: 200)
- `max_genes`: Maximum genes per cell (default: 5000)
- `min_cells`: Minimum cells per gene (default: 3)
- `max_mt_pct`: Maximum mitochondrial % (default: 20)

### Normalization
- `target_sum`: Normalization target (default: 10000)
- `n_hvg`: Number of HVGs (default: 2000)

### Clustering (in nextflow.config)
- `resolution`: Leiden resolution (default: [0.5, 1.0])
- `pca_dims`: PCA dimensions (default: 50)

### Marker Analysis
- `n_marker_genes`: Markers per cluster (default: 25)
- `marker_method`: Test method (default: wilcoxon)
- `n_plot_genes`: Genes in plots (default: 5)

## Usage

This is a hybrid pipeline with 3 sequential steps defined in `pipeline_order.json`.

### Step 1: QC + Normalization (Snakemake)
```bash
docker build -f Dockerfile.qc -t sc-qc .
docker run -v /path/to/input:/input:ro -v /path/to/output:/output sc-qc snakemake --cores 4
```

### Step 2: Integration + Clustering (Nextflow)
```bash
nextflow run nf-core/scdownstream -c nextflow.config --input /output/qc_normalized.h5ad
```

### Step 3: Marker Analysis (Snakemake)
```bash
docker build -f Dockerfile.marker -t sc-marker .
docker run -v /path/to/output:/input:ro -v /path/to/output:/output sc-marker snakemake --cores 4
```

## Requirements

- Docker
- Nextflow (for Step 2)
- Singularity (for nf-core containers)

## Author

PNU COLab