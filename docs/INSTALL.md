# Installation Guide

This guide provides detailed installation instructions for the WES Analysis Pipeline.

## System Requirements

### Minimum Requirements
- **Operating System**: Linux (Ubuntu 18.04+, CentOS 7+, RHEL 7+) or macOS 10.14+
- **Memory**: 16GB RAM (32GB+ recommended for large cohorts)
- **Storage**: 100GB free space minimum
  - Reference data: ~50GB
  - Example analysis: ~20GB per sample
  - Results vary by data size and number of samples
- **CPU**: 8+ cores recommended
- **Network**: Internet connection for downloading reference data and containers

### Recommended System Specifications
- **Memory**: 64GB+ RAM for optimal performance
- **Storage**: 500GB+ SSD storage
- **CPU**: 16+ cores with high clock speed
- **GPU**: Optional (for future GPU-accelerated variant calling)

## Software Dependencies

### Core Requirements

#### 1. Conda/Mamba (Required)
Install Miniconda or Anaconda:

```bash
# Download and install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Alternatively, install Mamba for faster package resolution
conda install -c conda-forge mamba
```

#### 2. Container Runtime (Required for DeepVariant)

**Option A: Docker (Recommended)**
```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install docker.io
sudo systemctl start docker
sudo systemctl enable docker

# Add user to docker group (logout/login required)
sudo usermod -aG docker $USER

# CentOS/RHEL
sudo yum install docker
sudo systemctl start docker
sudo systemctl enable docker
```

**Option B: Singularity (HPC environments)**
```bash
# Ubuntu/Debian
sudo apt-get install singularity-container

# CentOS/RHEL (requires EPEL)
sudo yum install epel-release
sudo yum install singularity
```

#### 3. Git (Usually pre-installed)
```bash
# Ubuntu/Debian
sudo apt-get install git

# CentOS/RHEL
sudo yum install git

# macOS
xcode-select --install
```

### Optional Dependencies

#### Java (for GATK - included in conda environment)
```bash
# If manual installation needed
sudo apt-get install openjdk-11-jdk
```

#### BWA, Samtools, etc. (installed via conda)
All bioinformatics tools are automatically installed via conda environment.

## Installation Steps

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/wes-analysis-pipeline.git
cd wes-analysis-pipeline
```

### 2. Create Conda Environment

```bash
# Create environment from file
conda env create -f environment.yml

# Or using mamba (faster)
mamba env create -f environment.yml

# Activate environment
conda activate wes-analysis
```

### 3. Verify Installation

```bash
# Test core tools
which bwa
which samtools
which gatk
which snakemake

# Check versions
bwa 2>&1 | head -3
samtools --version
gatk --version
snakemake --version
```

### 4. Test Container Runtime

**For Docker:**
```bash
# Test Docker installation
docker --version
docker run hello-world

# Test DeepVariant container pull
docker pull google/deepvariant:1.6.1
```

**For Singularity:**
```bash
# Test Singularity installation
singularity --version

# Test DeepVariant container
singularity pull docker://google/deepvariant:1.6.1
```

### 5. Prepare Reference Data

```bash
# Make script executable
chmod +x scripts/prepare_reference.sh

# Run reference preparation (requires ~1-2 hours and 50GB space)
./scripts/prepare_reference.sh
```

This will download and prepare:
- GRCh38 reference genome
- BWA, Samtools, and GATK indices
- Known variant sites (dbSNP, Mills indels, 1000G)
- SnpEff annotation database

## Installation Verification

### Quick Test

```bash
# Test with dry run (no execution)
./scripts/gatk_pipeline.sh -h
./scripts/deepvariant_pipeline.sh -h

# Test Snakemake workflows
./scripts/run_gatk_snakemake.sh -n
./scripts/run_deepvariant_snakemake.sh -n
```

### Full Validation

```bash
# Run installation test script
./tests/test_installation.sh

# Validate specific pipelines
./tests/test_gatk_pipeline.sh
./tests/test_deepvariant_pipeline.sh
```

## Troubleshooting

### Common Issues

#### 1. Conda Environment Creation Fails

**Issue**: Package conflicts or channel issues
```bash
# Solution: Update conda and try again
conda update conda
conda clean --all
conda env create -f environment.yml
```

#### 2. Docker Permission Denied

**Issue**: User not in docker group
```bash
# Solution: Add user to docker group
sudo usermod -aG docker $USER
# Logout and login again
```

#### 3. Container Runtime Not Working

**Issue**: Docker daemon not running
```bash
# Solution: Start Docker service
sudo systemctl start docker

# Check status
sudo systemctl status docker
```

#### 4. Reference Download Fails

**Issue**: Network connectivity or storage space
```bash
# Check available space
df -h

# Retry with specific reference
./scripts/prepare_reference.sh
```

#### 5. Memory Issues

**Issue**: Out of memory during execution
```bash
# Solution: Reduce thread count and increase memory allocation
./scripts/gatk_pipeline.sh -s sample1 -t 4 -m 32
```

### Getting Help

1. **Check logs**: All scripts generate detailed logs
2. **Review documentation**: See `docs/` directory
3. **Search issues**: Check GitHub issues
4. **Ask for help**: Create a new GitHub issue

## Performance Optimization

### System Configuration

#### 1. Increase File Descriptors
```bash
# Add to ~/.bashrc
ulimit -n 65536
```

#### 2. Optimize Filesystem
```bash
# Use faster filesystems (SSD preferred)
# Avoid network filesystems for temporary files
```

#### 3. Memory Settings
```bash
# Increase virtual memory for large datasets
echo 'vm.max_map_count=262144' | sudo tee -a /etc/sysctl.conf
sudo sysctl -p
```

### Container Optimization

#### Docker Settings
```bash
# Increase Docker memory allocation (in Docker Desktop)
# Recommended: 8GB+ for DeepVariant
```

#### Singularity Settings
```bash
# Set cache directory on fast storage
export SINGULARITY_CACHEDIR=/fast/storage/singularity_cache
```

## HPC Installation

### Module System

Many HPC systems use environment modules:

```bash
# Load required modules
module load conda/2023.03
module load singularity/3.8.0
module load java/11.0.2

# Then proceed with conda environment creation
```

### SLURM Integration

```bash
# Test SLURM integration
./scripts/run_gatk_snakemake.sh --slurm -j 100 -n
```

### Shared Storage

Ensure pipeline has access to:
- Shared conda environments
- Reference data on shared storage
- Results directory with appropriate permissions

## Next Steps

After successful installation:

1. **Read Usage Guide**: See [USAGE.md](USAGE.md)
2. **Configure Pipeline**: Edit `config/config.yaml`
3. **Prepare Input Data**: See [INPUT_FORMATS.md](INPUT_FORMATS.md)
4. **Run Test Analysis**: Use example data
5. **Scale to Production**: Configure for your HPC environment

## Support

For installation issues:
- Check [Troubleshooting Guide](TROUBLESHOOTING.md)
- Search [GitHub Issues](https://github.com/your-username/wes-analysis-pipeline/issues)
- Contact support team