#!/bin/bash

# validate_installation.sh
# Basic validation script for the WES Analysis Pipeline installation
# Used by CI/CD to verify that all components are working correctly

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_command() {
    local cmd="$1"
    local required="${2:-true}"

    if command -v "$cmd" &> /dev/null; then
        log_success "$cmd is available"
        return 0
    else
        if [ "$required" = "true" ]; then
            log_error "$cmd is not available and is required"
            return 1
        else
            log_warning "$cmd is not available but is optional"
            return 0
        fi
    fi
}

check_file() {
    local file="$1"
    local required="${2:-true}"

    if [ -f "$file" ]; then
        log_success "$file exists"
        return 0
    else
        if [ "$required" = "true" ]; then
            log_error "$file does not exist and is required"
            return 1
        else
            log_warning "$file does not exist but is optional"
            return 0
        fi
    fi
}

# Main validation
log_info "Starting WES Analysis Pipeline validation..."

# Check Python environment
log_info "Checking Python environment..."
check_command "python"
python --version

# Check core bioinformatics tools
log_info "Checking core bioinformatics tools..."
check_command "bwa"
check_command "samtools"
check_command "gatk"
check_command "fastqc"
check_command "multiqc"
check_command "snakemake"

# Check optional tools
log_info "Checking optional tools..."
check_command "bcftools" "false"
check_command "glnexus_cli" "false"
check_command "docker" "false"
check_command "singularity" "false"

# Check Python packages
log_info "Checking Python packages..."
python -c "
try:
    import pandas, numpy, matplotlib, seaborn, yaml
    print('✓ Core Python packages are available')
except ImportError as e:
    print(f'✗ Missing Python package: {e}')
    exit(1)
"

# Check configuration files
log_info "Checking configuration files..."
check_file "config/config.yaml"
check_file "config/cluster.yaml"
check_file "environment.yml"

# Validate YAML syntax
log_info "Validating YAML configuration syntax..."
python -c "
import yaml
try:
    with open('config/config.yaml') as f:
        yaml.safe_load(f)
    print('✓ config.yaml is valid')

    with open('config/cluster.yaml') as f:
        yaml.safe_load(f)
    print('✓ cluster.yaml is valid')

    with open('environment.yml') as f:
        yaml.safe_load(f)
    print('✓ environment.yml is valid')
except Exception as e:
    print(f'✗ YAML validation error: {e}')
    exit(1)
"

# Check pipeline scripts
log_info "Checking pipeline scripts..."
for script in scripts/*.sh; do
    if [ -f "$script" ]; then
        if [ -x "$script" ]; then
            log_success "$(basename "$script") is executable"
        else
            log_warning "$(basename "$script") is not executable"
        fi

        # Check script syntax
        if bash -n "$script"; then
            log_success "$(basename "$script") has valid syntax"
        else
            log_error "$(basename "$script") has syntax errors"
            exit 1
        fi
    fi
done

# Check Snakemake workflows
log_info "Checking Snakemake workflows..."
if [ -f "scripts/Snakefile_GATK" ]; then
    if snakemake -s scripts/Snakefile_GATK --lint 2>/dev/null; then
        log_success "GATK Snakefile syntax is valid"
    else
        log_warning "GATK Snakefile has lint warnings (may be acceptable)"
    fi
fi

if [ -f "scripts/Snakefile_DeepVariant" ]; then
    if snakemake -s scripts/Snakefile_DeepVariant --lint 2>/dev/null; then
        log_success "DeepVariant Snakefile syntax is valid"
    else
        log_warning "DeepVariant Snakefile has lint warnings (may be acceptable)"
    fi
fi

# Basic directory structure check
log_info "Checking directory structure..."
for dir in scripts config data reference results docs tests; do
    if [ -d "$dir" ]; then
        log_success "$dir/ directory exists"
    else
        log_warning "$dir/ directory does not exist"
    fi
done

# Summary
log_info "Validation completed!"
log_success "WES Analysis Pipeline installation appears to be working correctly"

echo
log_info "Next steps:"
echo "  1. Run './scripts/prepare_reference.sh' to download reference data"
echo "  2. Add your FASTQ files to the data/ directory"
echo "  3. Update config/config.yaml with your sample names"
echo "  4. Run the pipeline: './scripts/gatk_pipeline.sh -s your_sample'"

exit 0