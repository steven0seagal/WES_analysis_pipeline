#!/bin/bash

# deepvariant_pipeline.sh
# This script implements a WES pipeline using DeepVariant for variant calling
# From raw FASTQ files to annotated gVCF files
# Uses Docker/Singularity for DeepVariant execution

set -euo pipefail

# --- SCRIPT INFORMATION ---
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="DeepVariant WES Pipeline"

# --- CONFIGURATION ---
# Default paths - can be overridden by command line arguments
REF_DIR="reference"
REF_GENOME="${REF_DIR}/GRCh38.p13.genome.fa"
EXOME_TARGETS="${REF_DIR}/exome_targets.bed"

# Input files
DATA_DIR="data"
SAMPLE_NAME=""
FASTQ_R1=""
FASTQ_R2=""

# Output directory
RESULTS_DIR="results/deepvariant"

# Tool settings
THREADS=8
MEM_GB=16
DEEPVARIANT_VERSION="1.6.1"
CONTAINER_RUNTIME="docker"  # docker or singularity

# Read Group Information
RG_PL="ILLUMINA"
RG_CN="SequencingCenter"
RG_DS="WES"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# --- FUNCTIONS ---
usage() {
    cat << EOF
$SCRIPT_NAME v$SCRIPT_VERSION

Usage: $0 -s SAMPLE_NAME [OPTIONS]

Required Arguments:
  -s, --sample SAMPLE_NAME    Sample name (will look for SAMPLE_NAME_R1.fastq.gz and SAMPLE_NAME_R2.fastq.gz)

Optional Arguments:
  -d, --data-dir DIR          Input data directory (default: $DATA_DIR)
  -r, --ref-dir DIR           Reference directory (default: $REF_DIR)
  -o, --output-dir DIR        Output directory (default: $RESULTS_DIR)
  -t, --threads N             Number of threads (default: $THREADS)
  -m, --memory N              Memory in GB (default: $MEM_GB)
  -e, --exome-targets FILE    Exome targets BED file (default: $EXOME_TARGETS)
  --r1 FILE                   Specify R1 FASTQ file directly
  --r2 FILE                   Specify R2 FASTQ file directly
  --deepvariant-version VER   DeepVariant version (default: $DEEPVARIANT_VERSION)
  --container-runtime RT      Container runtime: docker or singularity (default: $CONTAINER_RUNTIME)
  --use-singularity           Use Singularity instead of Docker
  -h, --help                  Show this help message

Examples:
  $0 -s sample1
  $0 -s sample1 -t 16 -m 32
  $0 -s sample1 --use-singularity
  $0 -s sample1 --r1 /path/to/sample1_R1.fastq.gz --r2 /path/to/sample1_R2.fastq.gz

Notes:
  - This pipeline requires Docker or Singularity to run DeepVariant
  - BQSR is NOT performed as DeepVariant handles base quality errors internally
  - The pipeline is optimized for WES data using the WES model

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

run_fastqc() {
    local input_file="$1"
    local output_dir="$2"
    local log_file="$3"

    log_info "Running FastQC on $(basename "$input_file")"
    fastqc -o "$output_dir" -t "$THREADS" "$input_file" > "$log_file" 2>&1
}

check_container_runtime() {
    if [ "$CONTAINER_RUNTIME" = "docker" ]; then
        if ! command -v docker &> /dev/null; then
            log_error "Docker not found. Please install Docker or use --use-singularity option."
            exit 1
        fi
        # Check if Docker daemon is running
        if ! docker info &> /dev/null; then
            log_error "Docker daemon is not running. Please start Docker."
            exit 1
        fi
    elif [ "$CONTAINER_RUNTIME" = "singularity" ]; then
        if ! command -v singularity &> /dev/null; then
            log_error "Singularity not found. Please install Singularity or use Docker."
            exit 1
        fi
    else
        log_error "Invalid container runtime: $CONTAINER_RUNTIME. Use 'docker' or 'singularity'."
        exit 1
    fi
}

run_deepvariant_docker() {
    local input_bam="$1"
    local ref_genome="$2"
    local output_vcf="$3"
    local output_gvcf="$4"
    local targets_bed="$5"

    # Get absolute paths for Docker volume mounts
    local bam_abs=$(realpath "$input_bam")
    local ref_abs=$(realpath "$ref_genome")
    local targets_abs=$(realpath "$targets_bed")
    local output_dir_abs=$(realpath "$(dirname "$output_vcf")")

    # Check if BAM index exists
    if [ ! -f "${input_bam}.bai" ]; then
        log_error "BAM index not found: ${input_bam}.bai"
        exit 1
    fi

    log_info "Running DeepVariant via Docker..."
    docker run \
        --rm \
        -v "${bam_abs}:/input/$(basename "$input_bam"):ro" \
        -v "${bam_abs}.bai:/input/$(basename "$input_bam").bai:ro" \
        -v "${ref_abs}:/input/$(basename "$ref_genome"):ro" \
        -v "${ref_abs}.fai:/input/$(basename "$ref_genome").fai:ro" \
        -v "${targets_abs}:/input/$(basename "$targets_bed"):ro" \
        -v "${output_dir_abs}:/output" \
        "google/deepvariant:${DEEPVARIANT_VERSION}" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --ref="/input/$(basename "$ref_genome")" \
        --reads="/input/$(basename "$input_bam")" \
        --regions="/input/$(basename "$targets_bed")" \
        --output_vcf="/output/$(basename "$output_vcf")" \
        --output_gvcf="/output/$(basename "$output_gvcf")" \
        --num_shards="$THREADS" \
        --intermediate_results_dir="/output/intermediate_${SAMPLE_NAME}"
}

run_deepvariant_singularity() {
    local input_bam="$1"
    local ref_genome="$2"
    local output_vcf="$3"
    local output_gvcf="$4"
    local targets_bed="$5"

    # Download Singularity image if it doesn't exist
    local sif_file="deepvariant_${DEEPVARIANT_VERSION}.sif"
    if [ ! -f "$sif_file" ]; then
        log_info "Downloading DeepVariant Singularity image..."
        singularity pull "$sif_file" "docker://google/deepvariant:${DEEPVARIANT_VERSION}"
    fi

    log_info "Running DeepVariant via Singularity..."
    singularity exec \
        --bind "$(dirname "$input_bam"):/input" \
        --bind "$(dirname "$ref_genome"):/ref" \
        --bind "$(dirname "$targets_bed"):/targets" \
        --bind "$(dirname "$output_vcf"):/output" \
        "$sif_file" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --ref="/ref/$(basename "$ref_genome")" \
        --reads="/input/$(basename "$input_bam")" \
        --regions="/targets/$(basename "$targets_bed")" \
        --output_vcf="/output/$(basename "$output_vcf")" \
        --output_gvcf="/output/$(basename "$output_gvcf")" \
        --num_shards="$THREADS" \
        --intermediate_results_dir="/output/intermediate_${SAMPLE_NAME}"
}

# --- ARGUMENT PARSING ---
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -d|--data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        -r|--ref-dir)
            REF_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            RESULTS_DIR="$2"
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
        -e|--exome-targets)
            EXOME_TARGETS="$2"
            shift 2
            ;;
        --r1)
            FASTQ_R1="$2"
            shift 2
            ;;
        --r2)
            FASTQ_R2="$2"
            shift 2
            ;;
        --deepvariant-version)
            DEEPVARIANT_VERSION="$2"
            shift 2
            ;;
        --container-runtime)
            CONTAINER_RUNTIME="$2"
            shift 2
            ;;
        --use-singularity)
            CONTAINER_RUNTIME="singularity"
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
if [ -z "$SAMPLE_NAME" ]; then
    log_error "Sample name is required. Use -s or --sample option."
    usage
    exit 1
fi

# Set FASTQ files if not specified
if [ -z "$FASTQ_R1" ]; then
    FASTQ_R1="${DATA_DIR}/${SAMPLE_NAME}_R1.fastq.gz"
fi
if [ -z "$FASTQ_R2" ]; then
    FASTQ_R2="${DATA_DIR}/${SAMPLE_NAME}_R2.fastq.gz"
fi

# Update reference paths
REF_GENOME="${REF_DIR}/GRCh38.p13.genome.fa"

# Create output directories
SAMPLE_RESULTS_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_RESULTS_DIR}"/{bam,qc,vcf,gvcf,logs}
mkdir -p tmp

# Read Group Information
RG_ID="group1.${SAMPLE_NAME}"
RG_PU="unit1.${SAMPLE_NAME}"
RG_LB="lib1.${SAMPLE_NAME}"
RG_SM="${SAMPLE_NAME}"

# --- MAIN PIPELINE ---
log_info "Starting $SCRIPT_NAME v$SCRIPT_VERSION"
log_info "Sample: $SAMPLE_NAME"
log_info "Input R1: $FASTQ_R1"
log_info "Input R2: $FASTQ_R2"
log_info "Output directory: $SAMPLE_RESULTS_DIR"
log_info "Threads: $THREADS, Memory: ${MEM_GB}GB"
log_info "DeepVariant version: $DEEPVARIANT_VERSION"
log_info "Container runtime: $CONTAINER_RUNTIME"

# Check required tools
log_info "Checking required tools..."
for cmd in fastqc bwa samtools gatk; do
    check_command "$cmd"
done

# Check container runtime
check_container_runtime

# Check input files
log_info "Checking input files..."
check_file "$FASTQ_R1"
check_file "$FASTQ_R2"
check_file "$REF_GENOME"
check_file "$EXOME_TARGETS"

log_success "All checks passed. Starting pipeline..."

# --- STEP 1: Quality Control ---
log_info "=== STEP 1: Quality Control ==="
QC_DIR="${SAMPLE_RESULTS_DIR}/qc"
if [ ! -f "${QC_DIR}/${SAMPLE_NAME}_R1_fastqc.html" ] || [ ! -f "${QC_DIR}/${SAMPLE_NAME}_R2_fastqc.html" ]; then
    run_fastqc "$FASTQ_R1" "$QC_DIR" "${SAMPLE_RESULTS_DIR}/logs/fastqc_R1.log"
    run_fastqc "$FASTQ_R2" "$QC_DIR" "${SAMPLE_RESULTS_DIR}/logs/fastqc_R2.log"
    log_success "FastQC completed"
else
    log_info "FastQC already completed, skipping"
fi

# --- STEP 2: Alignment and Sorting ---
log_info "=== STEP 2: Alignment and Sorting ==="
ALIGNED_SORTED_BAM="${SAMPLE_RESULTS_DIR}/bam/${SAMPLE_NAME}.sorted.bam"

if [ ! -f "${ALIGNED_SORTED_BAM}" ]; then
    log_info "Aligning reads with BWA-MEM and sorting with Samtools..."
    bwa mem -t "${THREADS}" \
        -R "@RG\\tID:${RG_ID}\\tPL:${RG_PL}\\tPU:${RG_PU}\\tLB:${RG_LB}\\tSM:${RG_SM}\\tCN:${RG_CN}\\tDS:${RG_DS}" \
        "${REF_GENOME}" \
        "${FASTQ_R1}" \
        "${FASTQ_R2}" 2> "${SAMPLE_RESULTS_DIR}/logs/bwa_mem.log" | \
    samtools sort -@ "${THREADS}" -o "${ALIGNED_SORTED_BAM}" - 2> "${SAMPLE_RESULTS_DIR}/logs/samtools_sort.log"

    # Index the BAM file
    samtools index "${ALIGNED_SORTED_BAM}"
    log_success "Alignment and sorting completed"
else
    log_info "Aligned and sorted BAM already exists, skipping"
fi

# --- STEP 3: Mark Duplicates ---
log_info "=== STEP 3: Mark Duplicates ==="
MARKED_DUPS_BAM="${SAMPLE_RESULTS_DIR}/bam/${SAMPLE_NAME}.marked_dups.bam"
METRICS_FILE="${SAMPLE_RESULTS_DIR}/qc/${SAMPLE_NAME}.dup_metrics.txt"

if [ ! -f "${MARKED_DUPS_BAM}" ]; then
    log_info "Marking duplicates with Picard..."
    gatk --java-options "-Xmx${MEM_GB}G" MarkDuplicates \
        -I "${ALIGNED_SORTED_BAM}" \
        -O "${MARKED_DUPS_BAM}" \
        -M "${METRICS_FILE}" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        2> "${SAMPLE_RESULTS_DIR}/logs/mark_duplicates.log"
    log_success "Duplicate marking completed"
else
    log_info "Marked duplicates BAM already exists, skipping"
fi

# --- STEP 4: DeepVariant Variant Calling ---
log_info "=== STEP 4: DeepVariant Variant Calling ==="
log_warning "Note: BQSR is NOT performed for DeepVariant as it handles base quality errors internally"

OUTPUT_VCF="${SAMPLE_RESULTS_DIR}/vcf/${SAMPLE_NAME}.vcf.gz"
OUTPUT_GVCF="${SAMPLE_RESULTS_DIR}/gvcf/${SAMPLE_NAME}.g.vcf.gz"

if [ ! -f "${OUTPUT_GVCF}" ]; then
    log_info "Calling variants with DeepVariant..."

    if [ "$CONTAINER_RUNTIME" = "docker" ]; then
        run_deepvariant_docker "$MARKED_DUPS_BAM" "$REF_GENOME" "$OUTPUT_VCF" "$OUTPUT_GVCF" "$EXOME_TARGETS" \
            2> "${SAMPLE_RESULTS_DIR}/logs/deepvariant.log"
    elif [ "$CONTAINER_RUNTIME" = "singularity" ]; then
        run_deepvariant_singularity "$MARKED_DUPS_BAM" "$REF_GENOME" "$OUTPUT_VCF" "$OUTPUT_GVCF" "$EXOME_TARGETS" \
            2> "${SAMPLE_RESULTS_DIR}/logs/deepvariant.log"
    fi

    # Clean up intermediate files
    if [ -d "${SAMPLE_RESULTS_DIR}/vcf/intermediate_${SAMPLE_NAME}" ]; then
        rm -rf "${SAMPLE_RESULTS_DIR}/vcf/intermediate_${SAMPLE_NAME}"
    fi

    log_success "DeepVariant variant calling completed"
else
    log_info "DeepVariant gVCF already exists, skipping variant calling"
fi

# --- STEP 5: Quality Metrics ---
log_info "=== STEP 5: Quality Metrics ==="
HS_METRICS="${SAMPLE_RESULTS_DIR}/qc/${SAMPLE_NAME}.hs_metrics.txt"

if [ ! -f "${HS_METRICS}" ]; then
    log_info "Collecting hybrid selection metrics..."
    gatk --java-options "-Xmx${MEM_GB}G" CollectHsMetrics \
        -I "${MARKED_DUPS_BAM}" \
        -O "${HS_METRICS}" \
        -R "${REF_GENOME}" \
        -BAIT_INTERVALS "${EXOME_TARGETS}" \
        -TARGET_INTERVALS "${EXOME_TARGETS}" \
        2> "${SAMPLE_RESULTS_DIR}/logs/collect_hs_metrics.log"
    log_success "Quality metrics collection completed"
else
    log_info "Quality metrics already collected, skipping"
fi

# --- STEP 6: Basic Variant Statistics ---
log_info "=== STEP 6: Variant Statistics ==="
VARIANT_STATS="${SAMPLE_RESULTS_DIR}/qc/${SAMPLE_NAME}.variant_stats.txt"

if [ ! -f "${VARIANT_STATS}" ] && command -v bcftools &> /dev/null; then
    log_info "Generating variant statistics..."
    bcftools stats "${OUTPUT_VCF}" > "${VARIANT_STATS}" 2>/dev/null || true
    log_success "Variant statistics generated"
elif [ ! -f "${VARIANT_STATS}" ]; then
    log_warning "bcftools not found, skipping variant statistics"
fi

# --- STEP 7: Generate Summary Report ---
log_info "=== STEP 7: Generate Summary Report ==="
SUMMARY_REPORT="${SAMPLE_RESULTS_DIR}/${SAMPLE_NAME}_summary.txt"

cat > "${SUMMARY_REPORT}" << EOF
DeepVariant WES Pipeline Summary Report
=======================================
Sample: ${SAMPLE_NAME}
Pipeline Version: ${SCRIPT_VERSION}
Analysis Date: $(date)
Reference Genome: $(basename "${REF_GENOME}")
DeepVariant Version: ${DEEPVARIANT_VERSION}
Container Runtime: ${CONTAINER_RUNTIME}

Input Files:
- R1: ${FASTQ_R1}
- R2: ${FASTQ_R2}
- Exome Targets: ${EXOME_TARGETS}

Output Files:
- Preprocessed BAM: ${MARKED_DUPS_BAM}
- VCF: ${OUTPUT_VCF}
- gVCF: ${OUTPUT_GVCF}
- Quality Metrics: ${HS_METRICS}
- Duplicate Metrics: ${METRICS_FILE}

Processing Parameters:
- Threads: ${THREADS}
- Memory: ${MEM_GB}GB
- Model Type: WES (Whole Exome Sequencing)

Read Group Information:
- ID: ${RG_ID}
- Sample: ${RG_SM}
- Library: ${RG_LB}
- Platform: ${RG_PL}
- Platform Unit: ${RG_PU}

Key Differences from GATK Pipeline:
- NO Base Quality Score Recalibration (BQSR) performed
- DeepVariant handles base quality errors internally
- Uses deep learning CNN model for variant calling
- Simplified filtering requirements

Log Files:
- FastQC: ${SAMPLE_RESULTS_DIR}/logs/fastqc_*.log
- BWA-MEM: ${SAMPLE_RESULTS_DIR}/logs/bwa_mem.log
- Samtools: ${SAMPLE_RESULTS_DIR}/logs/samtools_sort.log
- Mark Duplicates: ${SAMPLE_RESULTS_DIR}/logs/mark_duplicates.log
- DeepVariant: ${SAMPLE_RESULTS_DIR}/logs/deepvariant.log
- Quality Metrics: ${SAMPLE_RESULTS_DIR}/logs/collect_hs_metrics.log

Next Steps:
1. Review quality metrics and FastQC reports
2. For cohort analysis, collect all gVCF files and run joint genotyping with GLnexus
3. Use scripts/glnexus_joint_calling.sh for cohort-level analysis

EOF

log_success "Summary report generated: ${SUMMARY_REPORT}"

# --- COMPLETION ---
log_success "=== PIPELINE COMPLETED SUCCESSFULLY ==="
log_info "Per-sample processing for ${SAMPLE_NAME} completed"
log_info "Total runtime: $SECONDS seconds"
log_info "Output gVCF: ${OUTPUT_GVCF}"
log_info "Output VCF: ${OUTPUT_VCF}"
log_info "Summary report: ${SUMMARY_REPORT}"
echo
log_info "To perform joint genotyping on multiple samples, use:"
log_info "  scripts/glnexus_joint_calling.sh -i results/deepvariant/*/gvcf/*.g.vcf.gz"

log_success "DeepVariant Pipeline finished successfully!"