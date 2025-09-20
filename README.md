# WES Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![Docker](https://img.shields.io/badge/docker-%230db7ed.svg)](https://www.docker.com/)
[![CI](https://github.com/yourusername/wes-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/yourusername/wes-pipeline/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)
[![GitHub issues](https://img.shields.io/github/issues/yourusername/wes-pipeline.svg)](https://github.com/yourusername/wes-pipeline/issues)
[![GitHub stars](https://img.shields.io/github/stars/yourusername/wes-pipeline.svg)](https://github.com/yourusername/wes-pipeline/stargazers)

A comprehensive Whole Exome Sequencing (WES) analysis pipeline repository with both GATK Best Practices and DeepVariant implementations. This repository provides production-ready pipelines for germline variant calling from raw FASTQ files to annotated VCFs.

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Pipeline Overview](#pipeline-overview)
- [Installation](#installation)
- [Usage](#usage)
- [Pipelines](#pipelines)
- [Configuration](#configuration)
- [Output](#output)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

## Features

- **Dual Pipeline Support**: Both GATK Best Practices and DeepVariant pipelines
- **Workflow Management**: Snakemake for reproducible and scalable execution
- **Container Support**: Docker/Singularity integration for DeepVariant
- **Quality Control**: Comprehensive QC with FastQC and MultiQC
- **Annotation**: SnpEff integration for variant annotation
- **Joint Calling**: Support for cohort-level analysis
- **CI/CD**: GitHub Actions for automated testing and deployment
- **Documentation**: Extensive guides and tutorials

## Quick Start

### Prerequisites

- Conda/Mamba for environment management
- Docker (for DeepVariant) or Singularity
- At least 16GB RAM, 8 CPU cores recommended

### Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/wes-pipeline.git
cd wes-pipeline
```

2. Set up the conda environment:
```bash
conda env create -f environment.yml
conda activate wes-analysis
```

3. Prepare reference data:
```bash
bash scripts/prepare_reference.sh
```

### Run a Pipeline

For GATK pipeline:
```bash
bash scripts/gatk_pipeline.sh -s sample1
```

For DeepVariant pipeline:
```bash
bash scripts/deepvariant_pipeline.sh -s sample1
```

## Pipeline Overview

### GATK Best Practices Pipeline

1. **Quality Control**: FastQC analysis of raw reads
2. **Alignment**: BWA-MEM alignment to reference genome
3. **Preprocessing**: Picard MarkDuplicates and Base Quality Score Recalibration (BQSR)
4. **Variant Calling**: GATK HaplotypeCaller in GVCF mode
5. **Annotation**: SnpEff functional annotation

### DeepVariant Pipeline

1. **Quality Control**: FastQC analysis of raw reads
2. **Alignment**: BWA-MEM alignment to reference genome
3. **Preprocessing**: Picard MarkDuplicates (no BQSR needed)
4. **Variant Calling**: DeepVariant via Docker/Singularity
5. **Annotation**: SnpEff functional annotation

## Installation

For detailed installation instructions, see [INSTALL.md](INSTALL.md).

## Usage

For detailed usage instructions, see [USAGE.md](USAGE.md).

## Pipelines

### Single Sample Analysis

Run individual pipelines on single samples:

```bash
# GATK pipeline
bash scripts/gatk_pipeline.sh -s sample1 -t 16 -m 32

# DeepVariant pipeline
bash scripts/deepvariant_pipeline.sh -s sample1 -t 16 -m 32
```

### Cohort Analysis

For multiple samples, use the joint calling scripts:

```bash
# GATK joint genotyping
bash scripts/joint_genotyping.sh -i results/gatk/*/gvcf/*.g.vcf.gz

# DeepVariant joint calling with GLnexus
bash scripts/glnexus_joint_calling.sh -i results/deepvariant/*/gvcf/*.g.vcf.gz
```

### Snakemake Workflows

Use Snakemake for automated parallel processing:

```bash
# GATK workflow
snakemake -s scripts/Snakefile_GATK --cores 16

# DeepVariant workflow (requires container setup)
snakemake -s scripts/Snakefile_DeepVariant --cores 16
```

## Configuration

All pipeline parameters are configured in `config/config.yaml`. Key settings include:

- Sample information and input files
- Reference genome and annotation files
- Computational resources (threads, memory)
- Tool-specific parameters
- Output options

See [docs/CONFIGURATION.md](docs/CONFIGURATION.md) for detailed configuration options.

## Output

### Directory Structure

```
results/
├── gatk/
│   └── sample1/
│       ├── bam/           # Aligned BAM files
│       ├── qc/            # Quality control reports
│       ├── recal/         # BQSR recalibration tables
│       ├── gvcf/          # gVCF files
│       └── logs/          # Pipeline logs
└── deepvariant/
    └── sample1/
        ├── bam/           # Aligned BAM files
        ├── qc/            # Quality control reports
        ├── gvcf/          # gVCF files
        └── logs/          # Pipeline logs
```

### Key Output Files

- **Analysis-ready BAM**: Recalibrated, duplicate-marked alignments
- **gVCF files**: Genomic VCFs for joint calling
- **Annotated VCF**: Functional annotations from SnpEff
- **Quality metrics**: Coverage, duplication, and variant statistics
- **MultiQC reports**: Aggregated quality control summaries

See [docs/OUTPUT_FORMATS.md](docs/OUTPUT_FORMATS.md) for detailed output descriptions.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Citation

If you use this pipeline in your research, please cite:

```
[DOI or citation information]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Documentation**: [docs/](docs/) directory
- **Issues**: [GitHub Issues](https://github.com/yourusername/wes-pipeline/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/wes-pipeline/discussions)

## Acknowledgments

This pipeline builds upon the excellent work of:

- [GATK](https://gatk.broadinstitute.org/) team
- [DeepVariant](https://github.com/google/deepvariant) developers
- [Snakemake](https://snakemake.github.io/) workflow management
- [Bioconda](https://bioconda.github.io/) for bioinformatics software distribution