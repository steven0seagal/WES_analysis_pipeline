#!/bin/bash

# joint_genotyping.sh
# Joint genotyping script for GATK gVCF files
# Performs cohort-level variant calling using GenomicsDBImport and GenotypeGVCFs

set -euo pipefail

# --- CONFIGURATION ---
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="GATK Joint Genotyping"

# Default parameters
REF_DIR="reference"
REF_GENOME="${REF_DIR}/GRCh38.p13.genome.fa"
EXOME_TARGETS="${REF_DIR}/exome_targets.bed"
OUTPUT_DIR="results/gatk/cohort"
COHORT_NAME="cohort"
THREADS=8
MEM_GB=32
GATK_JAVA_OPTS="-Xmx${MEM_GB}G -Djava.io.tmpdir=./tmp"

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
  -r, --ref-dir DIR           Reference directory (default: $REF_DIR)
  -e, --exome-targets FILE    Exome targets BED file (default: $EXOME_TARGETS)
  -t, --threads N             Number of threads (default: $THREADS)
  -m, --memory N              Memory in GB (default: $MEM_GB)
  -h, --help                  Show this help message

Examples:
  $0 -i "sample1.g.vcf.gz sample2.g.vcf.gz sample3.g.vcf.gz"
  $0 -i results/gatk/*/gvcf/*.g.vcf.gz -n my_cohort
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
        -r|--ref-dir)
            REF_DIR="$2"
            shift 2
            ;;
        -e|--exome-targets)
            EXOME_TARGETS="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEM_GB="$2"
            GATK_JAVA_OPTS="-Xmx${MEM_GB}G -Djava.io.tmpdir=./tmp"
            shift 2
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

# Update reference paths
REF_GENOME="${REF_DIR}/GRCh38.p13.genome.fa"

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p tmp

# Create a map file for GenomicsDBImport
SAMPLE_MAP="${OUTPUT_DIR}/${COHORT_NAME}_sample_map.txt"

log_info "Starting $SCRIPT_NAME v$SCRIPT_VERSION"
log_info "Cohort: $COHORT_NAME"
log_info "Number of samples: ${#GVCF_FILES[@]}"
log_info "Output directory: $OUTPUT_DIR"
log_info "Reference genome: $REF_GENOME"

# Check required tools
log_info "Checking required tools..."
check_command "gatk"

# Check reference files
log_info "Checking reference files..."
check_file "$REF_GENOME"
check_file "$EXOME_TARGETS"

# Validate and prepare gVCF files
log_info "Validating gVCF files and creating sample map..."
> "$SAMPLE_MAP"  # Clear the file

for i in "${!GVCF_FILES[@]}"; do
    gvcf_file="${GVCF_FILES[$i]}"

    if [ ! -f "$gvcf_file" ]; then
        log_error "gVCF file not found: $gvcf_file"
        exit 1
    fi

    # Extract sample name from filename
    sample_name=$(basename "$gvcf_file" .g.vcf.gz)

    # Add to sample map
    echo -e "${sample_name}\\t${gvcf_file}" >> "$SAMPLE_MAP"

    log_info "Sample $(($i + 1)): $sample_name -> $gvcf_file"
done

log_success "All gVCF files validated"

# --- STEP 1: GenomicsDBImport ---
log_info "=== STEP 1: GenomicsDBImport ==="
GENOMICS_DB="${OUTPUT_DIR}/${COHORT_NAME}_db"

if [ ! -d "$GENOMICS_DB" ]; then
    log_info "Creating GenomicsDB workspace..."

    # Build the -V arguments
    V_ARGS=""
    while IFS=$'\\t' read -r sample_name gvcf_path; do
        V_ARGS="$V_ARGS -V $gvcf_path"
    done < "$SAMPLE_MAP"

    gatk --java-options "${GATK_JAVA_OPTS}" GenomicsDBImport \
        $V_ARGS \
        --genomicsdb-workspace-path "$GENOMICS_DB" \
        -L "$EXOME_TARGETS" \
        --reader-threads 5 \
        --batch-size 50 \
        2> "${OUTPUT_DIR}/${COHORT_NAME}_genomicsdb_import.log"

    log_success "GenomicsDBImport completed"
else
    log_warning "GenomicsDB already exists, skipping import"
fi

# --- STEP 2: Joint Genotyping ---
log_info "=== STEP 2: Joint Genotyping ==="
COHORT_VCF="${OUTPUT_DIR}/${COHORT_NAME}.vcf.gz"

if [ ! -f "$COHORT_VCF" ]; then
    log_info "Performing joint genotyping..."

    gatk --java-options "${GATK_JAVA_OPTS}" GenotypeGVCFs \
        -R "$REF_GENOME" \
        -V "gendb://$GENOMICS_DB" \
        -O "$COHORT_VCF" \
        -L "$EXOME_TARGETS" \
        --include-non-variant-sites false \
        2> "${OUTPUT_DIR}/${COHORT_NAME}_genotype_gvcfs.log"

    log_success "Joint genotyping completed"
else
    log_warning "Cohort VCF already exists, skipping joint genotyping"
fi

# --- STEP 3: Basic Variant Statistics ---
log_info "=== STEP 3: Variant Statistics ==="
STATS_FILE="${OUTPUT_DIR}/${COHORT_NAME}_variant_stats.txt"

if [ ! -f "$STATS_FILE" ]; then
    log_info "Generating variant statistics..."

    gatk --java-options "${GATK_JAVA_OPTS}" VariantsToTable \
        -V "$COHORT_VCF" \
        -F CHROM -F POS -F TYPE -F QUAL -F AC -F AN -F AF \
        -O "${OUTPUT_DIR}/${COHORT_NAME}_variant_table.txt" \
        2> "${OUTPUT_DIR}/${COHORT_NAME}_variant_stats.log"

    # Generate summary statistics
    cat > "$STATS_FILE" << EOF
Joint Genotyping Summary Statistics
===================================
Cohort: ${COHORT_NAME}
Analysis Date: $(date)
Number of Samples: ${#GVCF_FILES[@]}

Output Files:
- Joint VCF: ${COHORT_VCF}
- Variant Table: ${OUTPUT_DIR}/${COHORT_NAME}_variant_table.txt
- GenomicsDB: ${GENOMICS_DB}

Sample Information:
$(cat "$SAMPLE_MAP")

Processing Parameters:
- Threads: ${THREADS}
- Memory: ${MEM_GB}GB
- Reference: ${REF_GENOME}
- Targets: ${EXOME_TARGETS}

Log Files:
- GenomicsDBImport: ${OUTPUT_DIR}/${COHORT_NAME}_genomicsdb_import.log
- GenotypeGVCFs: ${OUTPUT_DIR}/${COHORT_NAME}_genotype_gvcfs.log
- Variant Stats: ${OUTPUT_DIR}/${COHORT_NAME}_variant_stats.log

EOF

    log_success "Variant statistics generated"
else
    log_info "Variant statistics already exist, skipping"
fi

# --- STEP 4: Quality Summary ---
log_info "=== STEP 4: Quality Summary ==="

if command -v bcftools &> /dev/null; then
    log_info "Generating VCF summary with bcftools..."
    bcftools stats "$COHORT_VCF" > "${OUTPUT_DIR}/${COHORT_NAME}_bcftools_stats.txt" 2>/dev/null || true

    # Extract key metrics
    if [ -f "${OUTPUT_DIR}/${COHORT_NAME}_bcftools_stats.txt" ]; then
        echo >> "$STATS_FILE"
        echo "VCF Summary (bcftools stats):" >> "$STATS_FILE"
        echo "=============================" >> "$STATS_FILE"
        grep "^SN" "${OUTPUT_DIR}/${COHORT_NAME}_bcftools_stats.txt" | while read -r line; do
            echo "$line" | cut -f2- >> "$STATS_FILE"
        done
    fi
else
    log_warning "bcftools not found, skipping detailed VCF statistics"
fi

# --- COMPLETION ---
log_success "=== JOINT GENOTYPING COMPLETED SUCCESSFULLY ==="
log_info "Cohort analysis for ${COHORT_NAME} completed"
log_info "Number of samples processed: ${#GVCF_FILES[@]}"
log_info "Total runtime: $SECONDS seconds"
log_info "Final cohort VCF: ${COHORT_VCF}"
log_info "Summary statistics: ${STATS_FILE}"
echo
log_info "Next steps:"
log_info "1. Review variant statistics: ${STATS_FILE}"
log_info "2. Perform variant filtering and annotation"
log_info "3. Consider using GATK VQSR for large cohorts (>30 samples)"

log_success "Joint genotyping pipeline finished successfully!"