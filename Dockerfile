# WES Analysis Pipeline Docker Image
FROM continuumio/miniconda3:latest

# Set environment variables
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/conda/bin:$PATH

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools via conda
RUN conda install -c bioconda -c conda-forge \
    fastqc=0.12.1 \
    multiqc=1.24 \
    bwa=0.7.17 \
    samtools=1.19.2 \
    picard=3.1.1 \
    gatk4=4.5.0.0 \
    snpeff=5.2 \
    snakemake-minimal=8.14.0 \
    && conda clean --all -f -y

# Install Python dependencies
RUN pip install \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    pyyaml

# Create directories
RUN mkdir -p /opt/wes-pipeline /data /reference /results

# Copy pipeline files
COPY . /opt/wes-pipeline/
WORKDIR /opt/wes-pipeline

# Make scripts executable
RUN chmod +x scripts/*.sh

# Set up entrypoint
ENTRYPOINT ["/bin/bash"]
CMD ["--help"]

# Labels
LABEL maintainer="WES Pipeline Team" \
      description="Whole Exome Sequencing Analysis Pipeline" \
      version="1.0.0" \
      license="MIT"

# Default working directory for analysis
VOLUME ["/data", "/reference", "/results"]