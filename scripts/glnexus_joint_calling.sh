#!/bin/bash

# glnexus_joint_calling.sh
# Joint genotyping script for DeepVariant gVCF files using GLnexus
# Performs cohort-level variant calling optimized for DeepVariant output

set -euo pipefail

# --- CONFIGURATION ---
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="GLnexus Joint Calling for DeepVariant"

# Default parameters
OUTPUT_DIR="results/deepvariant/cohort"
COHORT_NAME="cohort"
THREADS=8
MEM_GB=32
CONFIG="DeepVariant_unfiltered"  # GLnexus configuration

# Input files
GVCF_FILES=()

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# --- FUNCTIONS ---
usage() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION

Usage: $0 -i GVCF_FILES [OPTIONS]

Required Arguments:
  -i, --input FILES           Space-separated list of gVCF files or directory with gVCF files

Optional Arguments:
  -o, --output-dir DIR        Output directory (default: $OUTPUT_DIR)
  -n, --cohort-name NAME      Cohort name prefix (default: $COHORT_NAME)
  -t, --threads N             Number of threads (default: $THREADS)
  -m, --memory N              Memory in GB (default: $MEM_GB)
  -c, --config CONFIG         GLnexus configuration (default: $CONFIG)
  --cleanup                   Remove intermediate GLnexus database after completion
  -h, --help                  Show this help message

GLnexus Configurations:
  - DeepVariant_unfiltered: Recommended for DeepVariant WES (default)
  - DeepVariantWGS: For whole genome sequencing data
  - gatk: For GATK gVCFs (use gatk joint genotyping instead)

Examples:
  $0 -i "sample1.g.vcf.gz sample2.g.vcf.gz sample3.g.vcf.gz"
  $0 -i results/deepvariant/*/gvcf/*.g.vcf.gz -n my_cohort
  $0 -i /path/to/gvcf/directory/ -t 16 -m 64

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

# --- ARGUMENT PARSING ---
CLEANUP=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            # Handle both directory and file list inputs
            input_arg="$2"
            if [ -d "$input_arg" ]; then
                # If it's a directory, find all gVCF files
                mapfile -t GVCF_FILES < <(find "$input_arg" -name "*.g.vcf.gz" -type f)
            else
                # If it's a file list, split by spaces
                IFS=' ' read -ra GVCF_FILES <<< "$input_arg"
            fi
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -n|--cohort-name)
            COHORT_NAME="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEM_GB="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        --cleanup)
            CLEANUP=true
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

# --- VALIDATION ---
if [ ${#GVCF_FILES[@]} -eq 0 ]; then
    log_error "No gVCF files specified. Use -i or --input option."
    usage
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

log_info "Starting $SCRIPT_NAME v$SCRIPT_VERSION"
log_info "Cohort: $COHORT_NAME"
log_info "Number of samples: ${#GVCF_FILES[@]}"
log_info "Output directory: $OUTPUT_DIR"
log_info "GLnexus configuration: $CONFIG"
log_info "Threads: $THREADS, Memory: ${MEM_GB}GB"

# Check required tools
log_info "Checking required tools..."
check_command "glnexus_cli"

# Validate gVCF files
log_info "Validating gVCF files..."
for i in "${!GVCF_FILES[@]}"; do
    gvcf_file="${GVCF_FILES[$i]}"

    if [ ! -f "$gvcf_file" ]; then
        log_error "gVCF file not found: $gvcf_file"
        exit 1
    fi

    # Extract sample name from filename
    sample_name=$(basename "$gvcf_file" .g.vcf.gz)
    log_info "Sample $(($i + 1)): $sample_name -> $gvcf_file"
done

log_success "All gVCF files validated"

# --- STEP 1: GLnexus Joint Calling ---
log_info "=== STEP 1: GLnexus Joint Calling ==="
COHORT_VCF="${OUTPUT_DIR}/${COHORT_NAME}.vcf.gz"
GLNEXUS_DB="${OUTPUT_DIR}/${COHORT_NAME}_glnexus.db"

if [ ! -f "$COHORT_VCF" ]; then
    log_info "Performing joint calling with GLnexus..."

    # Build the input arguments
    INPUT_ARGS=""
    for gvcf_file in "${GVCF_FILES[@]}"; do
        INPUT_ARGS="$INPUT_ARGS $gvcf_file"
    done

    # Run GLnexus
    glnexus_cli \
        --config "$CONFIG" \
        --threads "$THREADS" \
        --mem-gbytes "$MEM_GB" \
        --dir "$GLNEXUS_DB" \
        $INPUT_ARGS \
        2> "${OUTPUT_DIR}/${COHORT_NAME}_glnexus.log" | \
    bgzip -c > "$COHORT_VCF"

    # Index the VCF
    tabix -p vcf "$COHORT_VCF"

    log_success "GLnexus joint calling completed"
else
    log_warning "Cohort VCF already exists, skipping joint calling"
fi

# --- STEP 2: Basic Quality Filtering ---
log_info "=== STEP 2: Basic Quality Filtering ==="
FILTERED_VCF="${OUTPUT_DIR}/${COHORT_NAME}.filtered.vcf.gz"

if [ ! -f "$FILTERED_VCF" ] && command -v bcftools &> /dev/null; then
    log_info "Applying basic quality filters..."

    # Basic filtering for DeepVariant output
    # These are conservative filters - adjust based on your needs
    bcftools filter \
        -i 'QUAL>=20 && (INFO/DP>=10 && INFO/DP<=1000)' \
        -O z \
        -o "$FILTERED_VCF" \
        "$COHORT_VCF" \
        2> "${OUTPUT_DIR}/${COHORT_NAME}_filter.log"

    # Index the filtered VCF
    tabix -p vcf "$FILTERED_VCF"

    log_success "Basic quality filtering completed"
else
    if [ -f "$FILTERED_VCF" ]; then
        log_info "Filtered VCF already exists, skipping filtering"
    else
        log_warning "bcftools not found, skipping quality filtering"
    fi
fi

# --- STEP 3: Variant Statistics ---
log_info "=== STEP 3: Variant Statistics ==="
STATS_FILE="${OUTPUT_DIR}/${COHORT_NAME}_variant_stats.txt"

if [ ! -f "$STATS_FILE" ]; then
    log_info "Generating variant statistics..."

    cat > "$STATS_FILE" << EOF
GLnexus Joint Calling Summary Statistics
========================================
Cohort: ${COHORT_NAME}
Analysis Date: $(date)
Number of Samples: ${#GVCF_FILES[@]}
GLnexus Configuration: ${CONFIG}

Output Files:
- Joint VCF: ${COHORT_VCF}
- Filtered VCF: ${FILTERED_VCF}
- GLnexus Database: ${GLNEXUS_DB}

Sample Information:
EOF

    for i in "${!GVCF_FILES[@]}"; do
        sample_name=$(basename "${GVCF_FILES[$i]}" .g.vcf.gz)
        echo "  $((i + 1)). $sample_name: ${GVCF_FILES[$i]}" >> "$STATS_FILE"
    done

    cat >> "$STATS_FILE" << EOF

Processing Parameters:
- Threads: ${THREADS}
- Memory: ${MEM_GB}GB
- GLnexus Config: ${CONFIG}

Log Files:
- GLnexus: ${OUTPUT_DIR}/${COHORT_NAME}_glnexus.log
- Filtering: ${OUTPUT_DIR}/${COHORT_NAME}_filter.log

EOF

    log_success "Statistics file created"
else
    log_info "Statistics file already exists, skipping"
fi

# --- STEP 4: Quality Summary ---
log_info "=== STEP 4: Quality Summary ==="

if command -v bcftools &> /dev/null; then
    log_info "Generating VCF summary with bcftools..."

    # Stats for original VCF
    if [ -f "$COHORT_VCF" ]; then
        bcftools stats "$COHORT_VCF" > "${OUTPUT_DIR}/${COHORT_NAME}_bcftools_stats.txt" 2>/dev/null || true
    fi

    # Stats for filtered VCF
    if [ -f "$FILTERED_VCF" ]; then
        bcftools stats "$FILTERED_VCF" > "${OUTPUT_DIR}/${COHORT_NAME}_filtered_bcftools_stats.txt" 2>/dev/null || true
    fi

    # Extract key metrics
    if [ -f "${OUTPUT_DIR}/${COHORT_NAME}_bcftools_stats.txt" ]; then
        echo >> "$STATS_FILE"
        echo "Original VCF Summary (bcftools stats):" >> "$STATS_FILE"
        echo "=======================================" >> "$STATS_FILE"
        grep "^SN" "${OUTPUT_DIR}/${COHORT_NAME}_bcftools_stats.txt" | while read -r line; do
            echo "$line" | cut -f2- >> "$STATS_FILE"
        done
    fi

    if [ -f "${OUTPUT_DIR}/${COHORT_NAME}_filtered_bcftools_stats.txt" ]; then
        echo >> "$STATS_FILE"
        echo "Filtered VCF Summary (bcftools stats):" >> "$STATS_FILE"
        echo "=======================================" >> "$STATS_FILE"
        grep "^SN" "${OUTPUT_DIR}/${COHORT_NAME}_filtered_bcftools_stats.txt" | while read -r line; do
            echo "$line" | cut -f2- >> "$STATS_FILE"
        done
    fi
else
    log_warning "bcftools not found, skipping detailed VCF statistics"
fi

# --- STEP 5: Cleanup ---
if [ "$CLEANUP" = true ] && [ -d "$GLNEXUS_DB" ]; then
    log_info "=== STEP 5: Cleanup ==="
    log_info "Removing GLnexus database to save space..."
    rm -rf "$GLNEXUS_DB"
    log_success "GLnexus database removed"
fi

# --- COMPLETION ---
log_success "=== JOINT CALLING COMPLETED SUCCESSFULLY ==="
log_info "Cohort analysis for ${COHORT_NAME} completed"
log_info "Number of samples processed: ${#GVCF_FILES[@]}"
log_info "Total runtime: $SECONDS seconds"
log_info "Final cohort VCF: ${COHORT_VCF}"
if [ -f "$FILTERED_VCF" ]; then
    log_info "Filtered cohort VCF: ${FILTERED_VCF}"
fi
log_info "Summary statistics: ${STATS_FILE}"
echo
log_info "Next steps:"
log_info "1. Review variant statistics: ${STATS_FILE}"
log_info "2. Perform variant annotation with SnpEff or VEP"
log_info "3. Apply project-specific filtering criteria"
log_info "4. Consider DeepVariant-specific quality metrics for advanced filtering"

log_success "GLnexus joint calling pipeline finished successfully!"