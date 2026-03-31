FROM condaforge/miniforge3:latest

# Create environment with compatible Python version
RUN conda create -y -n scanpy python=3.11 && \
    conda clean -afy

# Activate environment and install packages
SHELL ["conda", "run", "-n", "scanpy", "/bin/bash", "-c"]

RUN conda install -y -c conda-forge -c bioconda \
    scanpy=1.10.2 \
    anndata=0.10.8 \
    snakemake-minimal \
    && conda clean -afy

# Set the environment to activate by default
ENV PATH=/opt/conda/envs/scanpy/bin:$PATH

# Setup pipeline
WORKDIR /pipeline
COPY Snakefile.qc Snakefile
COPY config.yaml .
COPY scripts/ scripts/

CMD ["snakemake", "--cores", "1"]