# Local Deployment Guide

This guide covers deploying the WES Analysis Pipeline on local machines and workstations.

## System Requirements

### Minimum Hardware

- **CPU**: 8 cores (16 recommended)
- **RAM**: 16 GB (32 GB recommended)
- **Storage**: 500 GB SSD
- **OS**: Linux (Ubuntu 18.04+, CentOS 7+)

### Recommended Hardware

- **CPU**: 16+ cores
- **RAM**: 64 GB
- **Storage**: 1+ TB NVMe SSD
- **Network**: 1 Gbps internet

## Installation Methods

### Method 1: Conda Environment (Recommended)

```bash
# Clone repository
git clone https://github.com/yourusername/wes-pipeline.git
cd wes-pipeline

# Create environment
conda env create -f environment.yml
conda activate wes-analysis

# Verify installation
bwa --version
gatk --version
```

### Method 2: Docker Container

```bash
# Build container
docker build -t wes-pipeline .

# Run container
docker run -it --rm \
  -v /path/to/data:/data \
  -v /path/to/reference:/reference \
  -v /path/to/results:/results \
  wes-pipeline

# Run pipeline inside container
bash scripts/gatk_pipeline.sh -s sample1 -d /data -o /results
```

### Method 3: Manual Installation

```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install openjdk-11-jdk python3-dev build-essential

# Install bioinformatics tools
conda install -c bioconda bwa samtools gatk4 picard snpeff
pip install snakemake pandas numpy

# Install FastQC
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip
sudo chmod +x FastQC/fastqc
sudo mv FastQC /usr/local/bin/
```

## Reference Data Setup

### Automated Setup

```bash
# Run reference preparation script
bash scripts/prepare_reference.sh
```

### Manual Setup

```bash
# Create reference directory
mkdir -p reference

# Download GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# Create indices
bwa index GCF_000001405.40_GRCh38.p14_genomic.fna
samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna
gatk CreateSequenceDictionary -R GCF_000001405.40_GRCh38.p14_genomic.fna
```

## Configuration

### Basic Configuration

```bash
# Edit config file
nano config/config.yaml

# Set paths
data_dir: "/path/to/your/data"
reference_genome: "/path/to/reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
exome_targets: "/path/to/exome_targets.bed"
```

### Resource Configuration

```yaml
# config/config.yaml
params:
  threads: 16
  mem_gb: 32

  # Adjust based on your system
  gatk:
    base_recalibrator_threads: 8
    haplotype_caller_threads: 8
```

## Running Pipelines

### Single Sample

```bash
# GATK pipeline
bash scripts/gatk_pipeline.sh -s sample1

# DeepVariant pipeline
bash scripts/deepvariant_pipeline.sh -s sample1
```

### Multiple Samples

```bash
# Process in parallel
parallel -j 4 'bash scripts/gatk_pipeline.sh -s {}' ::: sample1 sample2 sample3 sample4
```

### Snakemake Workflow

```bash
# Dry run
snakemake -s scripts/Snakefile_GATK --dry-run

# Execute
snakemake -s scripts/Snakefile_GATK --cores 16
```

## Monitoring and Troubleshooting

### Resource Monitoring

```bash
# Monitor CPU and memory
top -p $(pgrep -f gatk)

# Check disk usage
df -h

# Monitor I/O
iotop
```

### Common Issues

#### Out of Memory

```bash
# Reduce memory usage
export GATK_JAVA_OPTS="-Xmx16G"

# Use fewer threads
bash scripts/gatk_pipeline.sh -s sample1 -t 8 -m 16
```

#### Slow Performance

```bash
# Check CPU usage
mpstat 1

# Use SSD storage
# Ensure data is on fast storage
ls -la /path/to/data
```

#### Permission Issues

```bash
# Fix permissions
chmod -R 755 scripts/
chmod -R 755 results/
```

## Performance Optimization

### Hardware Optimization

1. **Use SSD Storage**: Place reference data and results on SSD
2. **Increase RAM**: More RAM allows larger thread counts
3. **CPU Cores**: Match thread count to available cores
4. **Network**: Fast internet for downloading reference data

### Software Optimization

```bash
# Optimize Java settings
export GATK_JAVA_OPTS="-Xmx32G -XX:+UseG1GC -XX:ParallelGCThreads=8"

# Use parallel processing
snakemake --cores all

# Enable compression
export BGZIP_OPT="-@ 8"
```

### Pipeline Optimization

```yaml
# config/config.yaml
output:
  keep_intermediate: false  # Save disk space
  compress_vcf: true        # Compress outputs

params:
  # Adjust based on data size
  bwa:
    threads: 16
  gatk:
    haplotype_caller_threads: 8
```

## Backup and Recovery

### Data Backup

```bash
# Backup important files
tar -czf backup_$(date +%Y%m%d).tar.gz \
  config/ \
  scripts/ \
  reference/ \
  results/

# Backup to external storage
rsync -av backup_*.tar.gz /external/drive/
```

### Recovery

```bash
# Restore from backup
tar -xzf backup_20231201.tar.gz

# Verify integrity
md5sum reference/*.fa
```

## Maintenance

### Regular Tasks

```bash
# Update conda environment
conda update --all

# Clean conda cache
conda clean --all

# Update reference data (if needed)
bash scripts/prepare_reference.sh

# Check disk space
df -h
du -sh results/
```

### Log Rotation

```bash
# Compress old logs
find results/ -name "*.log" -mtime +30 -exec gzip {} \;

# Clean old results (be careful!)
# find results/ -type f -mtime +90 -delete
```

## Scaling Up

### Multi-Core Processing

```bash
# Use all available cores
snakemake --cores all

# Limit concurrent jobs
snakemake --cores 32 --jobs 4
```

### Cluster Processing

For larger deployments, consider moving to HPC or cloud:

- **HPC**: Use SLURM/SGE job schedulers
- **Cloud**: AWS Batch, Google Cloud Life Sciences
- **Container Orchestration**: Kubernetes, Docker Swarm

See [HPC Deployment](DEPLOYMENT_HPC.md) and [Cloud Deployment](DEPLOYMENT_CLOUD.md) guides.