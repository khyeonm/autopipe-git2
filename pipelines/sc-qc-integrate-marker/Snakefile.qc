configfile: "config.yaml"

rule all:
    input:
        "/output/qc_normalized.h5ad"

rule qc_filter:
    """QC filtering: remove low-quality cells and genes."""
    input:
        h5ad=config["input_h5ad"]
    output:
        temp("/output/qc_filtered.h5ad")
    params:
        min_genes=config.get("min_genes", 200),
        min_cells=config.get("min_cells", 3),
        max_genes=config.get("max_genes", 5000),
        max_mt_pct=config.get("max_mt_pct", 20)
    log:
        "/output/logs/qc_filter.log"
    script:
        "scripts/qc_filter.py"

rule normalize:
    """Normalization and log transformation."""
    input:
        h5ad=rules.qc_filter.output
    output:
        "/output/qc_normalized.h5ad"
    params:
        target_sum=config.get("target_sum", 10000),
        n_top_genes=config.get("n_hvg", 2000)
    log:
        "/output/logs/normalize.log"
    script:
        "scripts/normalize.py"