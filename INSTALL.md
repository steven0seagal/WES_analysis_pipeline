# Installation Guide

This guide provides detailed instructions for installing and setting up the WES Analysis Pipeline.

## Table of Contents

- [System Requirements](#system-requirements)
- [Conda Environment Setup](#conda-environment-setup)
- [Docker Configuration](#docker-configuration)
- [Reference Data Preparation](#reference-data-preparation)
- [Pipeline Validation](#pipeline-validation)
- [Troubleshooting](#troubleshooting)

## System Requirements

### Minimum Requirements

- **Operating System**: Linux (Ubuntu 18.04+, CentOS 7+, or similar)
- **CPU**: 8 cores
- **RAM**: 16 GB
- **Storage**: 500 GB (for reference data and results)
- **Network**: Stable internet connection for downloading reference data

### Recommended Requirements

- **Operating System**: Linux (Ubuntu 20.04+ or CentOS 8+)
- **CPU**: 16+ cores
- **RAM**: 32+ GB
- **Storage**: 1+ TB SSD
- **Network**: High-speed internet (for large reference downloads)

### Software Prerequisites

- **Conda/Mamba**: For environment management
- **Git**: For version control
- **Docker** or **Singularity**: For DeepVariant execution
- **wget/curl**: For downloading reference data

## Conda Environment Setup

### Install Conda/Mamba

If you don't have Conda installed:

```bash
# Install Miniconda (recommended)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Or install Mamba (faster package resolution)
conda install mamba -n base -c conda-forge
```

### Create the Pipeline Environment

```bash
# Clone the repository
git clone https://github.com/yourusername/wes-pipeline.git
cd wes-pipeline

# Create the conda environment
conda env create -f environment.yml

# Activate the environment
conda activate wes-analysis

# Verify installation
which bwa samtools gatk fastqc
```

### Alternative: Manual Environment Creation

If you prefer to create the environment manually:

```bash
# Create environment
conda create -n wes-analysis python=3.11

# Activate
conda activate wes-analysis

# Install bioinformatics tools
conda install -c bioconda \
    fastqc=0.12.1 \
    multiqc=1.24 \
    bwa=0.7.17 \
    samtools=1.19.2 \
    picard=3.1.1 \
    gatk4=4.5.0.0 \
    snpeff=5.2 \
    snakemake-minimal=8.14.0

# Install additional tools
conda install -c conda-forge \
    pandas numpy matplotlib seaborn jupyter \
    black pre-commit \
    wget curl gzip
```

## Docker Configuration

### Install Docker

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install docker.io
sudo systemctl start docker
sudo systemctl enable docker
sudo usermod -aG docker $USER

# CentOS/RHEL
sudo yum install docker
sudo systemctl start docker
sudo systemctl enable docker
sudo usermod -aG docker $USER
```

### Test Docker Installation

```bash
# Test Docker
docker run hello-world

# Test DeepVariant container
docker run -it gcr.io/deepvariant-docker/deepvariant:1.6.1 /bin/bash
```

### Singularity (Alternative to Docker)

If Docker is not available (e.g., on HPC systems):

```bash
# Install Singularity
# Follow instructions at: https://sylabs.io/guides/3.7/user-guide/quick_start.html

# Pull DeepVariant container
singularity pull docker://gcr.io/deepvariant-docker/deepvariant:1.6.1
```

## Reference Data Preparation

### Automated Setup

Use the provided script to download and prepare reference data:

```bash
# Make script executable
chmod +x scripts/prepare_reference.sh

# Run the preparation script
bash scripts/prepare_reference.sh
```

This script will:
- Download GRCh38 reference genome
- Download known variant sites (dbSNP, Mills indels)
- Create BWA and Samtools indices
- Download SnpEff annotation databases

### Manual Setup

If you prefer manual setup or have existing reference data:

```bash
# Create reference directory
mkdir -p reference

# Download reference genome
cd reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38.p13.genome.fa

# Create indices
bwa index GRCh38.p13.genome.fa
samtools faidx GRCh38.p13.genome.fa
gatk CreateSequenceDictionary -R GRCh38.p13.genome.fa

# Download known variants
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz -O dbsnp_146.hg38.vcf.gz
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz.tbi -O dbsnp_146.hg38.vcf.gz.tbi

# Download Mills indels
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/structural_variant/GCF_000001405.40_GRCh38.p14_structural_variants.vcf.gz -O Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

# Download SnpEff database
snpeff download GRCh38.86
```

### Exome Targets

You must provide your own exome capture targets BED file:

```bash
# Place your exome targets file in the reference directory
cp /path/to/your/exome_targets.bed reference/
```

Common exome kits:
- Agilent SureSelect
- Illumina Nextera Rapid Capture
- IDT xGen

## Pipeline Validation

### Test Installation

Run the validation script:

```bash
bash scripts/validate_installation.sh
```

### Test with Example Data

```bash
# Download example data (if available)
# Run a small test pipeline
bash scripts/gatk_pipeline.sh -s test_sample --dry-run
```

### Check Tool Versions

```bash
# Verify all tools are accessible
bwa --version
samtools --version
gatk --version
fastqc --version
multiqc --version
snpeff -version
```

## Troubleshooting

### Common Issues

#### Conda Environment Issues

```bash
# Update conda
conda update conda

# Clean conda cache
conda clean --all

# Recreate environment
conda env remove -n wes-analysis
conda env create -f environment.yml
```

#### Docker Permission Issues

```bash
# Add user to docker group
sudo usermod -aG docker $USER

# Restart session or run:
newgrp docker
```

#### Reference Data Download Issues

```bash
# Check disk space
df -h

# Use alternative download method
wget --continue --tries=10 [URL]

# Or use curl
curl -C - -O [URL]
```

#### Memory Issues

```bash
# Check available memory
free -h

# Adjust pipeline parameters
# Edit config/config.yaml
# Reduce threads or memory usage
```

#### Permission Issues

```bash
# Make scripts executable
chmod +x scripts/*.sh

# Check file permissions
ls -la scripts/
```

### Getting Help

- Check the [troubleshooting guide](docs/TROUBLESHOOTING.md)
- Open an [issue](https://github.com/yourusername/wes-pipeline/issues)
- Review the [documentation](docs/)

### Performance Optimization

For large-scale analysis:

- Use SSD storage for reference data and results
- Increase RAM and CPU cores
- Use cluster/HPC for parallel processing
- Consider cloud deployment for elastic scaling

## Next Steps

After successful installation:

1. Review the [usage guide](USAGE.md)
2. Prepare your input data
3. Run test pipelines
4. Customize configuration as needed
5. Scale up for production analysis