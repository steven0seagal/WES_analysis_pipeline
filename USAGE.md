# Usage Guide

This guide provides detailed instructions for using the WES Analysis Pipeline.

## Table of Contents

- [Input Data Preparation](#input-data-preparation)
- [Configuration Options](#configuration-options)
- [Running Pipelines](#running-pipelines)
- [Output Interpretation](#output-interpretation)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)

## Input Data Preparation

### FASTQ Files

The pipeline expects paired-end FASTQ files with the following naming convention:

```
{sample_name}_R1.fastq.gz
{sample_name}_R2.fastq.gz
```

**Example:**
```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
```

### Directory Structure

Place your input files in the `data/` directory:

```
data/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

### Quality Requirements

- **Read Length**: 75-150 bp (typical for exome sequencing)
- **Quality Scores**: Q30+ for >80% of bases
- **Adapter Content**: <5% adapter contamination
- **Duplication Rate**: <20% (after PCR amplification)

### Exome Targets

You must provide an exome capture targets BED file:

```bash
# Place in reference directory
cp /path/to/exome_targets.bed reference/exome_targets.bed
```

The BED file should contain:
- Chromosome
- Start position (0-based)
- End position (1-based)
- Target name (optional)

**Example:**
```
chr1    69090   70008   target_1
chr1    69090   70008   target_2
```

## Configuration Options

### Main Configuration File

Edit `config/config.yaml` to customize pipeline behavior:

```yaml
# Sample information
samples:
  - sample1
  - sample2
  - sample3

# Directory paths
data_dir: "data"
ref_dir: "reference"
results_dir: "results"

# Computational resources
params:
  threads: 16
  mem_gb: 32

# Pipeline selection
run_gatk: true
run_deepvariant: true
```

### Tool-Specific Parameters

#### GATK Parameters

```yaml
gatk:
  base_recalibrator_threads: 8
  haplotype_caller_threads: 8
  joint_genotyping_threads: 16
```

#### DeepVariant Parameters

```yaml
deepvariant:
  version: "1.6.1"
  model_type: "WES"
  num_shards: 16
```

#### Quality Control

```yaml
fastqc:
  threads: 4
```

## Running Pipelines

### Single Sample Analysis

#### GATK Pipeline

```bash
# Basic usage
bash scripts/gatk_pipeline.sh -s sample1

# With custom resources
bash scripts/gatk_pipeline.sh -s sample1 -t 16 -m 32

# Specify input files directly
bash scripts/gatk_pipeline.sh -s sample1 --r1 /path/to/R1.fastq.gz --r2 /path/to/R2.fastq.gz

# Custom output directory
bash scripts/gatk_pipeline.sh -s sample1 -o results/custom_gatk
```

#### DeepVariant Pipeline

```bash
# Basic usage
bash scripts/deepvariant_pipeline.sh -s sample1

# With custom resources
bash scripts/deepvariant_pipeline.sh -s sample1 -t 16 -m 32

# Use Singularity instead of Docker
bash scripts/deepvariant_pipeline.sh -s sample1 --container-runtime singularity
```

### Multiple Samples

#### Batch Processing

```bash
# Process multiple samples
for sample in sample1 sample2 sample3; do
    bash scripts/gatk_pipeline.sh -s $sample
done
```

#### Snakemake Workflows

```bash
# GATK workflow
snakemake -s scripts/Snakefile_GATK --cores 16

# DeepVariant workflow
snakemake -s scripts/Snakefile_DeepVariant --cores 16

# Dry run to check workflow
snakemake -s scripts/Snakefile_GATK --dry-run

# Generate workflow diagram
snakemake -s scripts/Snakefile_GATK --dag | dot -Tpng > dag.png
```

### Cohort Analysis

#### GATK Joint Genotyping

```bash
# Collect all gVCF files
bash scripts/joint_genotyping.sh \
    -i results/gatk/*/gvcf/*.g.vcf.gz \
    -o results/gatk/cohort \
    -t 16
```

#### DeepVariant Joint Calling

```bash
# Use GLnexus for joint calling
bash scripts/glnexus_joint_calling.sh \
    -i results/deepvariant/*/gvcf/*.g.vcf.gz \
    -o results/deepvariant/cohort \
    -t 16
```

## Output Interpretation

### Directory Structure

```
results/
├── gatk/
│   └── sample1/
│       ├── bam/
│       │   ├── sample1.sorted.bam          # Aligned reads
│       │   ├── sample1.marked_dups.bam     # Duplicate-marked
│       │   └── sample1.analysis_ready.bam  # Final BAM
│       ├── qc/
│       │   ├── sample1_R1_fastqc.html      # FastQC reports
│       │   ├── sample1.dup_metrics.txt     # Duplication metrics
│       │   └── sample1.hs_metrics.txt      # Coverage metrics
│       ├── recal/
│       │   └── sample1.recal_data.table    # BQSR table
│       ├── gvcf/
│       │   └── sample1.g.vcf.gz            # Genomic VCF
│       └── logs/                           # Pipeline logs
└── deepvariant/
    └── sample1/
        ├── bam/
        │   ├── sample1.sorted.bam          # Aligned reads
        │   ├── sample1.marked_dups.bam     # Duplicate-marked
        │   └── sample1.analysis_ready.bam  # Final BAM
        ├── qc/
        │   ├── sample1_R1_fastqc.html      # FastQC reports
        │   └── sample1.dup_metrics.txt     # Duplication metrics
        ├── gvcf/
        │   └── sample1.g.vcf.gz            # Genomic VCF
        └── logs/                           # Pipeline logs
```

### Key Metrics to Review

#### Coverage Metrics

Check `hs_metrics.txt` for:
- **Mean Coverage**: Should be 50-100x for exome
- **Coverage Uniformity**: >80% at 10x coverage
- **Target Coverage**: >95% of targets covered at 20x

#### Quality Metrics

- **Duplication Rate**: <20% acceptable
- **Base Quality**: Q30+ for >80% of bases
- **Mapping Quality**: >95% reads with MQ>20

#### Variant Calling

- **Ti/Tv Ratio**: ~2.0-2.1 for WES
- **Het/Hom Ratio**: ~1.5-2.0
- **Novel Variants**: Check against known databases

### MultiQC Reports

Generate comprehensive QC reports:

```bash
# Install MultiQC if not in environment
pip install multiqc

# Generate report
multiqc results/ -o results/multiqc_report/
```

## Advanced Usage

### Custom Reference Genome

```bash
# Use different reference
echo "reference_genome: 'reference/custom.fa'" >> config/config.yaml

# Update indices
bwa index reference/custom.fa
samtools faidx reference/custom.fa
gatk CreateSequenceDictionary -R reference/custom.fa
```

### Cluster Execution

#### SLURM

```yaml
# config/cluster.yaml
__default__:
  partition: normal
  time: "24:00:00"
  mem: "16G"
  cpus: 8
```

```bash
# Submit to cluster
snakemake -s scripts/Snakefile_GATK \
    --cluster-config config/cluster.yaml \
    --cluster "sbatch -p {cluster.partition} -t {cluster.time} -c {cluster.cpus} --mem={cluster.mem}" \
    --jobs 10
```

#### SGE

```yaml
# config/cluster.yaml
__default__:
  queue: all.q
  pe: smp
  time: "24:00:00"
```

### Pipeline Customization

#### Skip Steps

```bash
# Skip BQSR (not recommended for GATK)
# Edit pipeline script or Snakemake rule

# Skip certain QC steps
# Comment out rules in Snakefile
```

#### Custom Annotation

```bash
# Add custom annotation databases
snpeff build -gtf22 -v GRCh38.86 /path/to/custom.gtf

# Use different annotation source
# Edit snpeff configuration in pipeline
```

### Performance Optimization

#### Memory Management

```yaml
# config/config.yaml
params:
  java_opts: "-Xmx64G -Djava.io.tmpdir=./tmp"
  threads: 32
```

#### Parallel Processing

```bash
# Use all available cores
snakemake -s scripts/Snakefile_GATK --cores all

# Limit concurrent jobs
snakemake -s scripts/Snakefile_GATK --cores 16 --jobs 4
```

## Troubleshooting

### Common Issues

#### Pipeline Fails at Alignment

```bash
# Check FASTQ format
zcat sample_R1.fastq.gz | head -20

# Verify reference indices
ls -la reference/*.fai reference/*.dict reference/*.sa
```

#### Memory Errors

```bash
# Reduce Java heap size
export GATK_JAVA_OPTS="-Xmx16G"

# Increase system memory
# Or reduce threads
```

#### Docker Issues

```bash
# Check Docker service
sudo systemctl status docker

# Pull latest DeepVariant image
docker pull gcr.io/deepvariant-docker/deepvariant:1.6.1
```

#### File Permission Issues

```bash
# Fix permissions
chmod -R 755 scripts/
chmod -R 755 results/
```

### Getting Help

- Check pipeline logs in `results/*/logs/`
- Review FastQC and MultiQC reports
- Validate input files format
- Check tool versions and dependencies

### Performance Monitoring

```bash
# Monitor resource usage
top -p $(pgrep -f gatk)

# Check disk I/O
iotop

# Monitor memory usage
free -h
```

## Next Steps

- Review quality metrics and optimize parameters
- Perform cohort analysis with joint calling
- Annotate variants and filter results
- Integrate with downstream analysis tools
- Set up automated pipelines for production use