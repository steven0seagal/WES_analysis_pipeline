#!/bin/bash

# test_environment_locally.sh
# Local testing script to validate environment before pushing to CI/CD
# This helps catch issues before they cause CI failures

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

# Test environment creation
test_environment_creation() {
    log_info "Testing conda environment creation..."

    # Check if environment already exists
    if conda env list | grep -q "wes-analysis-test"; then
        log_warning "Test environment already exists, removing it..."
        conda env remove -n wes-analysis-test -y
    fi

    # Try to create environment
    log_info "Creating test environment from environment.yml..."
    if conda env create -f environment.yml -n wes-analysis-test; then
        log_success "Environment created successfully"

        # Test activation
        log_info "Testing environment activation..."
        if conda activate wes-analysis-test 2>/dev/null; then
            log_success "Environment activated successfully"

            # Test core tools
            log_info "Testing core tools availability..."

            # Test each tool
            local tools=("python" "bwa" "samtools" "gatk" "fastqc" "multiqc" "snakemake")
            local failed_tools=()

            for tool in "${tools[@]}"; do
                if command -v "$tool" &> /dev/null; then
                    log_success "$tool is available"
                else
                    log_error "$tool is not available"
                    failed_tools+=("$tool")
                fi
            done

            # Test Python packages
            log_info "Testing Python packages..."
            if python -c "import pandas, numpy, matplotlib, seaborn, yaml" 2>/dev/null; then
                log_success "Core Python packages are available"
            else
                log_error "Some Python packages are missing"
                failed_tools+=("python_packages")
            fi

            # Clean up test environment
            conda deactivate 2>/dev/null || true
            log_info "Removing test environment..."
            conda env remove -n wes-analysis-test -y

            if [ ${#failed_tools[@]} -eq 0 ]; then
                log_success "All tools are working correctly"
                return 0
            else
                log_error "Failed tools: ${failed_tools[*]}"
                return 1
            fi
        else
            log_error "Failed to activate environment"
            return 1
        fi
    else
        log_error "Failed to create environment"
        return 1
    fi
}

# Test YAML syntax
test_yaml_syntax() {
    log_info "Testing YAML syntax..."

    local yaml_files=(
        "environment.yml"
        "config/config.yaml"
        "config/cluster.yaml"
        ".pre-commit-config.yaml"
    )

    for yaml_file in "${yaml_files[@]}"; do
        if [ -f "$yaml_file" ]; then
            if python -c "import yaml; yaml.safe_load(open('$yaml_file'))" 2>/dev/null; then
                log_success "$yaml_file syntax is valid"
            else
                log_error "$yaml_file has syntax errors"
                return 1
            fi
        else
            log_warning "$yaml_file does not exist"
        fi
    done

    return 0
}

# Test script syntax
test_script_syntax() {
    log_info "Testing shell script syntax..."

    local scripts_dir="scripts"
    local failed_scripts=()

    if [ -d "$scripts_dir" ]; then
        for script in "$scripts_dir"/*.sh; do
            if [ -f "$script" ]; then
                if bash -n "$script"; then
                    log_success "$(basename "$script") syntax is valid"
                else
                    log_error "$(basename "$script") has syntax errors"
                    failed_scripts+=("$(basename "$script")")
                fi
            fi
        done

        if [ ${#failed_scripts[@]} -eq 0 ]; then
            log_success "All shell scripts have valid syntax"
            return 0
        else
            log_error "Scripts with syntax errors: ${failed_scripts[*]}"
            return 1
        fi
    else
        log_error "Scripts directory does not exist"
        return 1
    fi
}

# Test Snakemake workflow syntax
test_snakemake_syntax() {
    log_info "Testing Snakemake workflow syntax..."

    local workflows=(
        "scripts/Snakefile_GATK"
        "scripts/Snakefile_DeepVariant"
    )

    # Need to activate environment for snakemake
    if conda env list | grep -q "wes-analysis"; then
        conda activate wes-analysis 2>/dev/null || true
    fi

    for workflow in "${workflows[@]}"; do
        if [ -f "$workflow" ]; then
            if snakemake -s "$workflow" --lint 2>/dev/null; then
                log_success "$(basename "$workflow") syntax is valid"
            else
                log_warning "$(basename "$workflow") has lint warnings (may be acceptable)"
            fi
        else
            log_warning "$workflow does not exist"
        fi
    done

    conda deactivate 2>/dev/null || true
    return 0
}

# Test directory structure
test_directory_structure() {
    log_info "Testing directory structure..."

    local required_dirs=(
        "scripts"
        "config"
        "docs"
        ".github/workflows"
    )

    local optional_dirs=(
        "data"
        "reference"
        "results"
        "tests"
        "notebooks"
        "tutorials"
    )

    for dir in "${required_dirs[@]}"; do
        if [ -d "$dir" ]; then
            log_success "$dir/ directory exists"
        else
            log_error "$dir/ directory is missing"
            return 1
        fi
    done

    for dir in "${optional_dirs[@]}"; do
        if [ -d "$dir" ]; then
            log_success "$dir/ directory exists"
        else
            log_info "$dir/ directory does not exist (optional)"
        fi
    done

    return 0
}

# Main function
main() {
    log_info "Starting local environment testing..."
    echo "=================================================="

    local failed_tests=()

    # Run tests
    echo
    if ! test_directory_structure; then
        failed_tests+=("directory_structure")
    fi

    echo
    if ! test_yaml_syntax; then
        failed_tests+=("yaml_syntax")
    fi

    echo
    if ! test_script_syntax; then
        failed_tests+=("script_syntax")
    fi

    echo
    if ! test_snakemake_syntax; then
        failed_tests+=("snakemake_syntax")
    fi

    echo
    if ! test_environment_creation; then
        failed_tests+=("environment_creation")
    fi

    # Summary
    echo
    echo "=================================================="
    log_info "Local testing summary:"

    if [ ${#failed_tests[@]} -eq 0 ]; then
        log_success "🎉 All tests passed! Environment should work in CI/CD"
        echo
        log_info "Next steps:"
        echo "  1. Commit and push your changes"
        echo "  2. Check CI/CD pipeline results"
        echo "  3. Fix any remaining issues"
        exit 0
    else
        log_error "❌ Some tests failed: ${failed_tests[*]}"
        echo
        log_info "Please fix the issues above before pushing to CI/CD"
        exit 1
    fi
}

# Check if conda is available
if ! command -v conda &> /dev/null; then
    log_error "Conda is not available. Please install conda/miniconda first."
    exit 1
fi

# Run main function
main "$@"