#!/bin/bash

# run_gatk_snakemake.sh
# Wrapper script to run the GATK Snakemake workflow
# Provides easy command-line interface and validation

set -euo pipefail

# --- Configuration ---
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="GATK Snakemake Workflow Runner"

# Default parameters
CONFIG_FILE="config/config.yaml"
SNAKEFILE="scripts/Snakefile_GATK"
CORES=8
DRYRUN=false
UNLOCK=false
FORCE_RERUN=false
REPORT=false

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# --- Functions ---
usage() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION

Usage: $0 [OPTIONS]

Options:
  -c, --config FILE           Configuration file (default: $CONFIG_FILE)
  -s, --snakefile FILE        Snakefile (default: $SNAKEFILE)
  -j, --cores N               Number of cores (default: $CORES)
  -n, --dry-run               Perform a dry run
  -u, --unlock                Unlock working directory
  -f, --force                 Force rerun of all rules
  -r, --report                Generate HTML report
  --cluster CMD               Submit jobs to cluster using CMD
  --slurm                     Use SLURM cluster profile
  -h, --help                  Show this help message

Examples:
  $0                          # Run with default settings
  $0 -j 16 -n                # Dry run with 16 cores
  $0 --slurm -j 100          # Submit to SLURM cluster
  $0 -r                      # Generate HTML report

EOF
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $1"
}

check_file() {
    if [ ! -f "$1" ]; then
        log_error "File not found: $1"
        exit 1
    fi
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "Command '$1' not found. Please ensure it's installed and in PATH."
        exit 1
    fi
}

# --- Argument Parsing ---
CLUSTER_CMD=""
USE_SLURM=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -s|--snakefile)
            SNAKEFILE="$2"
            shift 2
            ;;
        -j|--cores)
            CORES="$2"
            shift 2
            ;;
        -n|--dry-run)
            DRYRUN=true
            shift
            ;;
        -u|--unlock)
            UNLOCK=true
            shift
            ;;
        -f|--force)
            FORCE_RERUN=true
            shift
            ;;
        -r|--report)
            REPORT=true
            shift
            ;;
        --cluster)
            CLUSTER_CMD="$2"
            shift 2
            ;;
        --slurm)
            USE_SLURM=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# --- Validation ---
log_info "Starting $SCRIPT_NAME v$SCRIPT_VERSION"

# Check required files
check_file "$CONFIG_FILE"
check_file "$SNAKEFILE"

# Check required commands
check_command "snakemake"

if [ "$USE_SLURM" = true ]; then
    check_command "sbatch"
    check_command "squeue"
fi

# --- Main Execution ---

# Unlock if requested
if [ "$UNLOCK" = true ]; then
    log_info "Unlocking Snakemake working directory..."
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" --unlock
    log_success "Directory unlocked"
    exit 0
fi

# Generate report if requested
if [ "$REPORT" = true ]; then
    log_info "Generating Snakemake HTML report..."
    snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" --report results/gatk_workflow_report.html
    log_success "Report generated: results/gatk_workflow_report.html"
    exit 0
fi

# Build Snakemake command
SNAKE_CMD="snakemake"
SNAKE_CMD+=" --snakefile $SNAKEFILE"
SNAKE_CMD+=" --configfile $CONFIG_FILE"
SNAKE_CMD+=" --cores $CORES"
SNAKE_CMD+=" --use-conda"
SNAKE_CMD+=" --conda-frontend mamba"
SNAKE_CMD+=" --printshellcmds"
SNAKE_CMD+=" --reason"

# Add optional parameters
if [ "$DRYRUN" = true ]; then
    SNAKE_CMD+=" --dry-run"
    log_info "Performing dry run..."
fi

if [ "$FORCE_RERUN" = true ]; then
    SNAKE_CMD+=" --forceall"
    log_warning "Forcing rerun of all rules"
fi

# Cluster configuration
if [ "$USE_SLURM" = true ]; then
    log_info "Configuring for SLURM cluster..."
    SNAKE_CMD+=" --cluster 'sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={threads} --job-name={rule}.{wildcards}'"
    SNAKE_CMD+=" --cluster-config config/cluster.yaml"
    SNAKE_CMD+=" --jobs 100"
elif [ -n "$CLUSTER_CMD" ]; then
    log_info "Using custom cluster command: $CLUSTER_CMD"
    SNAKE_CMD+=" --cluster '$CLUSTER_CMD'"
    SNAKE_CMD+=" --jobs 100"
fi

# Display configuration
log_info "Configuration:"
echo "  Config file: $CONFIG_FILE"
echo "  Snakefile: $SNAKEFILE"
echo "  Cores: $CORES"
echo "  Dry run: $DRYRUN"
echo "  Force rerun: $FORCE_RERUN"
echo "  Use SLURM: $USE_SLURM"

if [ -n "$CLUSTER_CMD" ]; then
    echo "  Cluster command: $CLUSTER_CMD"
fi

echo

# Validate workflow first
log_info "Validating workflow..."
if snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" --dry-run --quiet; then
    log_success "Workflow validation passed"
else
    log_error "Workflow validation failed"
    exit 1
fi

# Display workflow summary
log_info "Workflow summary:"
snakemake --snakefile "$SNAKEFILE" --configfile "$CONFIG_FILE" --summary

echo

# Execute workflow
log_info "Executing Snakemake workflow..."
log_info "Command: $SNAKE_CMD"

if eval "$SNAKE_CMD"; then
    log_success "Workflow completed successfully!"

    if [ "$DRYRUN" = false ]; then
        echo
        log_info "Output files:"
        echo "  - Final annotated VCF: results/gatk/cohort/cohort.annotated.vcf.gz"
        echo "  - MultiQC report: results/gatk/multiqc/multiqc_report.html"
        echo "  - Individual sample summaries: results/gatk/*/sample_summary.txt"
        echo
        log_info "Next steps:"
        echo "  1. Review the MultiQC report for quality assessment"
        echo "  2. Examine individual sample summaries"
        echo "  3. Perform variant filtering and downstream analysis"
    fi
else
    log_error "Workflow execution failed"
    exit 1
fi