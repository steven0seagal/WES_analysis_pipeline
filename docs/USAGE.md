# Usage Guide

This guide provides comprehensive instructions for running the WES Analysis Pipeline.

## Overview

The pipeline offers multiple execution modes:
- **Bash Scripts**: Simple, educational, single-sample execution
- **Snakemake Workflows**: Production-ready, parallel, multi-sample execution
- **Joint Calling**: Cohort-level variant calling for multiple samples

## Prerequisites

Before running any pipeline:

1. **Complete Installation**: Follow [INSTALL.md](INSTALL.md)
2. **Prepare Reference Data**: Run `./scripts/prepare_reference.sh`
3. **Activate Environment**: `conda activate wes-analysis`
4. **Prepare Input Data**: See [Input Data Preparation](#input-data-preparation)

## Input Data Preparation

### Required Input Files

#### 1. FASTQ Files
```bash
# Expected naming convention:
data/
├── sample1_R1.fastq.gz    # Forward reads
├── sample1_R2.fastq.gz    # Reverse reads
├── sample2_R1.fastq.gz
├── sample2_R2.fastq.gz
└── ...
```

#### 2. Exome Target Regions (BED file)
```bash
# User must provide exome capture kit coordinates
reference/exome_targets.bed
```

Example BED format:
```
chr1    65419    65433
chr1    65520    65573
chr1    69037    70008
```

#### 3. Sample Configuration
Edit `config/config.yaml`:
```yaml
samples:
  - sample1
  - sample2
  - sample3
```

## Execution Modes

### 1. Bash Script Execution

#### GATK Pipeline

**Basic Usage:**
```bash
./scripts/gatk_pipeline.sh -s sample1
```

**Advanced Usage:**
```bash
./scripts/gatk_pipeline.sh \
    -s sample1 \
    -t 16 \
    -m 32 \
    -d /path/to/data \
    -r /path/to/reference \
    -o /path/to/output
```

**Parameters:**
- `-s, --sample`: Sample name (required)
- `-t, --threads`: Number of threads (default: 8)
- `-m, --memory`: Memory in GB (default: 16)
- `-d, --data-dir`: Input data directory
- `-r, --ref-dir`: Reference directory
- `-o, --output-dir`: Output directory
- `--r1, --r2`: Direct FASTQ file paths

#### DeepVariant Pipeline

**Basic Usage:**
```bash
./scripts/deepvariant_pipeline.sh -s sample1
```

**With Singularity:**
```bash
./scripts/deepvariant_pipeline.sh -s sample1 --use-singularity
```

**Advanced Usage:**
```bash
./scripts/deepvariant_pipeline.sh \
    -s sample1 \
    -t 16 \
    -m 32 \
    --deepvariant-version 1.6.1 \
    --container-runtime docker
```

**Additional Parameters:**
- `--deepvariant-version`: DeepVariant version
- `--container-runtime`: docker or singularity
- `--use-singularity`: Use Singularity instead of Docker

### 2. Snakemake Workflow Execution

#### GATK Workflow

**Basic Usage:**
```bash
./scripts/run_gatk_snakemake.sh -j 16
```

**Production Usage:**
```bash
./scripts/run_gatk_snakemake.sh \
    -c config/config.yaml \
    -j 32 \
    --cluster-config config/cluster.yaml
```

**Dry Run (Recommended First):**
```bash
./scripts/run_gatk_snakemake.sh -n
```

#### DeepVariant Workflow

**Basic Usage:**
```bash
./scripts/run_deepvariant_snakemake.sh -j 16
```

**With Singularity:**
```bash
./scripts/run_deepvariant_snakemake.sh -j 16 --use-singularity
```

**Cluster Execution:**
```bash
./scripts/run_deepvariant_snakemake.sh --slurm -j 100
```

### 3. Joint Calling (Cohort Analysis)

#### GATK Joint Genotyping

After running GATK pipeline on all samples:

```bash
./scripts/joint_genotyping.sh \
    -i "results/gatk/*/gvcf/*.g.vcf.gz" \
    -o results/gatk/cohort \
    -n my_cohort \
    -t 16 \
    -m 64
```

#### DeepVariant Joint Calling (GLnexus)

After running DeepVariant pipeline on all samples:

```bash
./scripts/glnexus_joint_calling.sh \
    -i "results/deepvariant/*/gvcf/*.g.vcf.gz" \
    -o results/deepvariant/cohort \
    -n my_cohort \
    -t 16 \
    -m 64
```

## Cluster Execution

### SLURM Integration

#### Automatic SLURM Submission
```bash
# GATK workflow on SLURM
./scripts/run_gatk_snakemake.sh --slurm -j 100

# DeepVariant workflow on SLURM
./scripts/run_deepvariant_snakemake.sh --slurm -j 100
```

#### Custom SLURM Configuration
Edit `config/cluster.yaml`:
```yaml
__default__:
  partition: "normal"
  time: "02:00:00"
  mem: "8G"

deepvariant:
  time: "06:00:00"
  mem: "16G"
  partition: "gpu"
```

#### Manual SLURM Job Submission
```bash
# Submit as SLURM job
sbatch --job-name=wes_gatk \
       --partition=normal \
       --time=24:00:00 \
       --mem=64G \
       --cpus-per-task=16 \
       scripts/run_gatk_snakemake.sh -j 16
```

### SGE Integration

```bash
./scripts/run_gatk_snakemake.sh \
    --cluster "qsub -q all.q -pe smp {threads} -l mem_req={resources.mem_mb}M"
```

## Configuration Options

### Main Configuration (`config/config.yaml`)

#### Sample Information
```yaml
samples:
  - sample1
  - sample2
  - sample3

# OR load from file
samples: !include samples.txt
```

#### Resource Allocation
```yaml
params:
  threads: 16
  mem_gb: 32
  java_opts: "-Xmx32G -Djava.io.tmpdir=./tmp"
```

#### Tool Parameters
```yaml
params:
  bwa:
    threads: 16
  gatk:
    base_recalibrator_threads: 8
    haplotype_caller_threads: 8
  deepvariant:
    version: "1.6.1"
    model_type: "WES"
    num_shards: 16
```

#### Output Options
```yaml
output:
  keep_intermediate: false
  compress_vcf: true
  generate_multiqc: true
```

## Monitoring and Debugging

### Log Files

All scripts generate comprehensive logs:

```bash
# GATK pipeline logs
results/gatk/{sample}/logs/
├── bwa_mem.log
├── mark_duplicates.log
├── base_recalibrator.log
├── apply_bqsr.log
└── haplotype_caller.log

# DeepVariant pipeline logs
results/deepvariant/{sample}/logs/
├── bwa_mem.log
├── mark_duplicates.log
└── deepvariant.log
```

### Snakemake Monitoring

```bash
# Check workflow status
snakemake -s scripts/Snakefile_GATK --summary

# Generate workflow report
snakemake -s scripts/Snakefile_GATK --report results/workflow_report.html

# Monitor cluster jobs
snakemake -s scripts/Snakefile_GATK --cluster-status scripts/cluster_status.py
```

### Quality Control

#### MultiQC Reports
```bash
# Automatically generated in:
results/{pipeline}/multiqc/multiqc_report.html
```

#### Manual QC Checks
```bash
# Check sample summaries
cat results/gatk/{sample}/{sample}_summary.txt
cat results/deepvariant/{sample}/{sample}_summary.txt

# Check variant statistics
bcftools stats results/gatk/cohort/cohort.vcf.gz
bcftools stats results/deepvariant/cohort/cohort.vcf.gz
```

## Performance Optimization

### Resource Allocation

#### Memory Recommendations
- **Single sample**: 16-32GB
- **Joint calling**: 32-64GB
- **Large cohorts**: 64-128GB

#### Thread Allocation
- **BWA alignment**: 8-16 threads
- **GATK HaplotypeCaller**: 4-8 threads
- **DeepVariant**: 8-32 threads (matches num_shards)

### Storage Optimization

#### Temporary Files
```bash
# Use fast local storage for temporary files
export TMPDIR=/fast/local/tmp
mkdir -p $TMPDIR
```

#### Intermediate File Cleanup
```yaml
# In config.yaml
output:
  keep_intermediate: false  # Delete intermediate BAMs
  cleanup_temp: true        # Clean temporary files
```

### Container Optimization

#### Docker Settings
```bash
# Allocate sufficient memory to Docker
# Docker Desktop: Preferences > Resources > Memory (8GB+)

# Use local Docker registry for faster pulls
docker pull google/deepvariant:1.6.1
```

#### Singularity Cache
```bash
# Set Singularity cache on fast storage
export SINGULARITY_CACHEDIR=/fast/storage/singularity
```

## Troubleshooting Common Issues

### Pipeline Failures

#### Out of Memory
```bash
# Reduce threads and increase memory per thread
./scripts/gatk_pipeline.sh -s sample1 -t 4 -m 32
```

#### Container Issues
```bash
# Test container manually
docker run google/deepvariant:1.6.1 /opt/deepvariant/bin/run_deepvariant --help

# Check container logs
docker logs <container_id>
```

#### File Permission Errors
```bash
# Fix permissions
chmod -R 755 results/
chmod +x scripts/*.sh
```

### Snakemake Issues

#### Workflow Corruption
```bash
# Unlock and restart
./scripts/run_gatk_snakemake.sh --unlock
./scripts/run_gatk_snakemake.sh -j 16
```

#### Missing Files
```bash
# Force regeneration
./scripts/run_gatk_snakemake.sh -j 16 --force
```

#### Cluster Connection Issues
```bash
# Test cluster connectivity
sinfo  # For SLURM
qstat  # For SGE
```

## Best Practices

### 1. Start Small
- Test with one sample first
- Use dry runs to validate workflows
- Check intermediate outputs

### 2. Monitor Resources
- Use system monitoring tools (`htop`, `iostat`)
- Check disk space regularly
- Monitor memory usage

### 3. Version Control
- Document software versions
- Keep configuration files in version control
- Track analysis parameters

### 4. Quality Control
- Always review MultiQC reports
- Check sample summaries
- Validate variant call statistics

### 5. Reproducibility
- Use identical software versions
- Document all parameters
- Save configuration files with results

## Example Workflows

### Single Sample Analysis

```bash
# 1. Activate environment
conda activate wes-analysis

# 2. Run GATK pipeline
./scripts/gatk_pipeline.sh -s sample1 -t 16 -m 32

# 3. Run DeepVariant pipeline
./scripts/deepvariant_pipeline.sh -s sample1 -t 16 -m 32

# 4. Compare results
./scripts/compare_pipelines.sh sample1
```

### Cohort Analysis

```bash
# 1. Configure samples
nano config/config.yaml  # Add all sample names

# 2. Run Snakemake workflow
./scripts/run_gatk_snakemake.sh -j 32

# 3. Joint calling
./scripts/joint_genotyping.sh -i "results/gatk/*/gvcf/*.g.vcf.gz"

# 4. Review results
firefox results/gatk/multiqc/multiqc_report.html
```

### HPC Production Run

```bash
# 1. Configure for cluster
cp config/cluster.yaml.template config/cluster.yaml
nano config/cluster.yaml

# 2. Submit to cluster
./scripts/run_gatk_snakemake.sh --slurm -j 200

# 3. Monitor progress
squeue -u $USER
watch -n 30 'squeue -u $USER'

# 4. Check results when complete
ls results/gatk/cohort/
```

## Next Steps

After successful execution:

1. **Review Output**: See [OUTPUT_FORMATS.md](OUTPUT_FORMATS.md)
2. **Variant Filtering**: Apply quality filters
3. **Annotation**: Add functional annotations
4. **Downstream Analysis**: Population genetics, association studies
5. **Visualization**: Generate plots and reports

## Support

For usage questions:
- Check [FAQ](FAQ.md)
- Review [Troubleshooting Guide](TROUBLESHOOTING.md)
- Search [GitHub Issues](https://github.com/your-username/wes-analysis-pipeline/issues)
- Post in [Discussions](https://github.com/your-username/wes-analysis-pipeline/discussions)