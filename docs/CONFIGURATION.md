# Configuration Guide

This guide explains all configuration options available in the WES Analysis Pipeline.

## Configuration Files Overview

The pipeline uses several configuration files:

- **`config/config.yaml`**: Main configuration file
- **`config/cluster.yaml`**: HPC cluster settings
- **`environment.yml`**: Software dependencies
- **`requirements.txt`**: Python dependencies

## Main Configuration (`config/config.yaml`)

### Sample Configuration

#### Basic Sample List
```yaml
# List of sample names
samples:
  - sample1
  - sample2
  - sample3
```

#### Sample List from File
```yaml
# Load samples from external file
samples: !include samples.txt
```

Where `samples.txt` contains:
```
sample1
sample2
sample3
```

#### Sample Metadata
```yaml
# Advanced sample configuration with metadata
samples:
  sample1:
    condition: "case"
    batch: "batch1"
    sex: "male"
  sample2:
    condition: "control"
    batch: "batch1"
    sex: "female"
```

### Directory Paths

```yaml
# Input and output directories
data_dir: "data"                    # Input FASTQ directory
ref_dir: "reference"                # Reference data directory
results_dir: "results"              # Output directory
scripts_dir: "scripts"              # Pipeline scripts

# Specific file paths
reference_genome: "reference/GRCh38.p13.genome.fa"
dbsnp: "reference/dbsnp_146.hg38.vcf.gz"
mills_indels: "reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
exome_targets: "reference/exome_targets.bed"
```

### Computational Resources

#### Basic Resource Allocation
```yaml
params:
  threads: 8          # Default thread count
  mem_gb: 16          # Default memory in GB
  java_opts: "-Xmx16G -Djava.io.tmpdir=./tmp"
```

#### Tool-Specific Resource Settings
```yaml
params:
  # FastQC settings
  fastqc:
    threads: 2

  # BWA alignment settings
  bwa:
    threads: 16
    mem_sort_threads: 8

  # GATK-specific settings
  gatk:
    base_recalibrator_threads: 8
    haplotype_caller_threads: 8
    joint_genotyping_threads: 16

  # DeepVariant settings
  deepvariant:
    version: "1.6.1"
    model_type: "WES"           # WES, WGS, or PACBIO
    num_shards: 16              # Parallel processing shards
```

### Read Group Information

```yaml
params:
  rg:
    platform: "ILLUMINA"         # Sequencing platform
    library_prefix: "lib"        # Library naming prefix
    platform_unit_prefix: "unit" # Platform unit prefix
```

This generates read groups like:
```
@RG\tID:unit.sample1\tSM:sample1\tPL:ILLUMINA\tLB:libsample1
```

### Pipeline Selection

```yaml
# Enable/disable specific pipelines
run_gatk: true
run_deepvariant: true

# Cohort name for joint calling
cohort_name: "my_cohort"
```

### Output Options

```yaml
output:
  # Keep intermediate files (useful for debugging)
  keep_intermediate: false

  # Compress output VCF files
  compress_vcf: true

  # Generate MultiQC report
  generate_multiqc: true

  # Clean up temporary files
  cleanup_temp: true
```

### Advanced Container Options

```yaml
advanced:
  # Use container for DeepVariant (required)
  use_deepvariant_container: true

  # Container runtime: docker or singularity
  container_runtime: "docker"

  # Temporary directory
  tmp_dir: "./tmp"

  # Clean up temporary files automatically
  cleanup_temp: true
```

### Annotation Settings

```yaml
params:
  snpeff:
    genome_db: "GRCh38.86"       # SnpEff database version
    memory: "8g"                 # Memory allocation for SnpEff
```

## Cluster Configuration (`config/cluster.yaml`)

### SLURM Configuration

#### Default Settings
```yaml
__default__:
  partition: "normal"                    # Default partition
  time: "02:00:00"                      # Default time limit
  mem: "8G"                             # Default memory
  output: "logs/slurm/{rule}_{wildcards}_%j.out"
  error: "logs/slurm/{rule}_{wildcards}_%j.err"
```

#### Rule-Specific Settings
```yaml
# Memory-intensive steps
mark_duplicates:
  time: "02:00:00"
  mem: "16G"

base_recalibrator:
  time: "04:00:00"
  mem: "16G"

haplotype_caller:
  time: "08:00:00"
  mem: "32G"

# Long-running steps
deepvariant:
  time: "12:00:00"
  mem: "16G"
  partition: "gpu"              # Use GPU partition if available

genomicsdb_import:
  time: "06:00:00"
  mem: "64G"

genotype_gvcfs:
  time: "12:00:00"
  mem: "64G"

# Quick steps
fastqc_raw:
  time: "00:30:00"
  mem: "4G"

multiqc:
  time: "00:30:00"
  mem: "4G"
```

### SGE Configuration

```yaml
__default__:
  queue: "all.q"
  pe: "smp"
  h_vmem: "8G"
  h_rt: "02:00:00"

deepvariant:
  h_vmem: "16G"
  h_rt: "12:00:00"
  queue: "gpu.q"
```

## Performance Profiles

### Profile for Small Datasets (< 10 samples)
```yaml
# config/performance_profiles/small.yaml
params:
  threads: 8
  mem_gb: 16
  bwa:
    threads: 8
  gatk:
    haplotype_caller_threads: 4
  deepvariant:
    num_shards: 8
```

### Profile for Medium Datasets (10-50 samples)
```yaml
# config/performance_profiles/medium.yaml
params:
  threads: 16
  mem_gb: 32
  bwa:
    threads: 16
  gatk:
    haplotype_caller_threads: 8
    joint_genotyping_threads: 16
  deepvariant:
    num_shards: 16
```

### Profile for Large Datasets (50+ samples)
```yaml
# config/performance_profiles/large.yaml
params:
  threads: 32
  mem_gb: 64
  bwa:
    threads: 32
  gatk:
    haplotype_caller_threads: 16
    joint_genotyping_threads: 32
  deepvariant:
    num_shards: 32
```

## Environment Configuration

### Conda Environment (`environment.yml`)

#### Core Dependencies
```yaml
name: wes-analysis
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  # Quality control
  - fastqc=0.12.1
  - multiqc=1.24

  # Alignment and processing
  - bwa=0.7.17
  - samtools=1.19.2
  - picard=3.1.1

  # Variant calling
  - gatk4=4.5.0.0

  # Annotation
  - snpeff=5.2
  - bcftools=1.19

  # Workflow management
  - snakemake-minimal=8.14.0

  # Joint calling for DeepVariant
  - glnexus=1.4.1
```

#### Development Dependencies
```yaml
  # Development tools
  - black=23.12.0
  - pre-commit=3.6.0
  - pytest=7.4.0
```

### Python Requirements (`requirements.txt`)

```txt
# Core scientific computing
pandas>=2.1.0
numpy>=1.24.0
matplotlib>=3.8.0
seaborn>=0.13.0

# Workflow management
snakemake>=8.0.0
pyyaml>=6.0

# Development tools
black>=23.0.0
pre-commit>=3.6.0

# Documentation
sphinx>=7.0.0
sphinx-rtd-theme>=1.3.0
```

## Configuration Templates

### Minimal Configuration
```yaml
# Minimal config for testing
samples:
  - test_sample

params:
  threads: 4
  mem_gb: 8

run_gatk: true
run_deepvariant: false
```

### Production Configuration
```yaml
# Production config with all options
samples: !include samples.txt

data_dir: "/data/fastq"
ref_dir: "/reference/GRCh38"
results_dir: "/results/wes_analysis"

params:
  threads: 32
  mem_gb: 64
  java_opts: "-Xmx64G -Djava.io.tmpdir=/tmp"

  bwa:
    threads: 32
    mem_sort_threads: 16

  gatk:
    base_recalibrator_threads: 16
    haplotype_caller_threads: 16
    joint_genotyping_threads: 32

  deepvariant:
    version: "1.6.1"
    model_type: "WES"
    num_shards: 32

output:
  keep_intermediate: false
  compress_vcf: true
  generate_multiqc: true
  cleanup_temp: true

advanced:
  container_runtime: "singularity"
  tmp_dir: "/tmp/wes_pipeline"

run_gatk: true
run_deepvariant: true
cohort_name: "production_cohort"
```

### HPC Configuration
```yaml
# HPC-optimized configuration
samples: !include samples.txt

# Use shared storage paths
data_dir: "/shared/data/fastq"
ref_dir: "/shared/reference/GRCh38"
results_dir: "/shared/results/wes_analysis"

params:
  threads: 16
  mem_gb: 32

advanced:
  container_runtime: "singularity"
  tmp_dir: "/tmp/$USER/wes_pipeline"
  use_shared_conda: true

# Cluster-specific settings
cluster:
  slurm:
    partition: "normal"
    time: "24:00:00"
    mem_per_cpu: "4G"
    account: "your_account"
```

## Configuration Validation

### Schema Validation
```yaml
# config/schema.yaml (optional validation schema)
type: object
properties:
  samples:
    type: array
    items:
      type: string
    minItems: 1
  params:
    type: object
    properties:
      threads:
        type: integer
        minimum: 1
        maximum: 256
      mem_gb:
        type: integer
        minimum: 4
        maximum: 1024
```

### Configuration Testing

```bash
# Test configuration validity
./scripts/validate_config.py config/config.yaml

# Dry run to check configuration
./scripts/run_gatk_snakemake.sh -n
```

## Environment Variables

### Runtime Environment Variables
```bash
# Set temporary directory
export TMPDIR=/fast/tmp

# Set Singularity cache
export SINGULARITY_CACHEDIR=/shared/singularity_cache

# Set conda environment location
export CONDA_ENVS_PATH=/shared/conda/envs

# Set Java options globally
export JAVA_OPTS="-Xmx32G"
```

### Pipeline-Specific Variables
```bash
# Override configuration paths
export WES_CONFIG_DIR=/custom/config
export WES_REFERENCE_DIR=/custom/reference
export WES_RESULTS_DIR=/custom/results

# Container settings
export DEEPVARIANT_VERSION=1.6.1
export CONTAINER_RUNTIME=singularity
```

## Best Practices

### 1. Resource Planning
- Start with conservative resource allocation
- Monitor actual usage and adjust
- Consider node specifications on HPC systems

### 2. Path Management
- Use absolute paths for shared systems
- Keep configuration files under version control
- Document any custom paths or settings

### 3. Scalability
- Test with small datasets first
- Scale resources based on data size
- Consider parallel execution limits

### 4. Reproducibility
- Pin software versions in configuration
- Document all custom settings
- Keep configuration with analysis results

### 5. Security
- Avoid hardcoding passwords or tokens
- Use environment variables for sensitive data
- Restrict file permissions appropriately

## Troubleshooting Configuration Issues

### Common Problems

#### Invalid YAML Syntax
```bash
# Check YAML syntax
python -c "import yaml; yaml.safe_load(open('config/config.yaml'))"
```

#### Resource Allocation Errors
```bash
# Check available resources
free -h          # Memory
nproc            # CPU cores
df -h            # Disk space
```

#### Path Resolution Issues
```bash
# Test path accessibility
ls -la config/config.yaml
ls -la reference/GRCh38.p13.genome.fa
```

#### Container Configuration Issues
```bash
# Test container runtime
docker --version
singularity --version

# Test container access
docker run hello-world
singularity run docker://hello-world
```

### Configuration Debugging

```bash
# Enable debug mode
export WES_DEBUG=1

# Verbose configuration checking
./scripts/run_gatk_snakemake.sh -n -v

# Print effective configuration
snakemake -s scripts/Snakefile_GATK --configfile config/config.yaml --print-compilation
```

## Support

For configuration help:
- Review [Installation Guide](INSTALL.md)
- Check [Usage Guide](USAGE.md)
- See [Troubleshooting Guide](TROUBLESHOOTING.md)
- Ask in [GitHub Discussions](https://github.com/your-username/wes-analysis-pipeline/discussions)