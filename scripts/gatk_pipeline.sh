#!/bin/bash

# gatk_pipeline.sh
# This script implements the GATK Best Practices workflow for Whole Exome Sequencing
# From raw FASTQ files to analysis-ready BAM and gVCF files
# Includes quality control, alignment, preprocessing, and variant calling

set -euo pipefail  # Exit on error, unset variable, or pipe failure

# --- SCRIPT INFORMATION ---
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="GATK WES Pipeline"

# --- CONFIGURATION ---
# Default paths - can be overridden by command line arguments
REF_DIR="reference"
REF_GENOME="${REF_DIR}/GRCh38.p13.genome.fa"
DBSNP_VCF="${REF_DIR}/dbsnp_146.hg38.vcf.gz"
MILLS_INDELS_VCF="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
THOUSAND_GENOMES_VCF="${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
EXOME_TARGETS="${REF_DIR}/exome_targets.bed"

# Input files
DATA_DIR="data"
SAMPLE_NAME=""
FASTQ_R1=""
FASTQ_R2=""

# Output directory
RESULTS_DIR="results/gatk"

# Tool settings
THREADS=8
MEM_GB=16
GATK_JAVA_OPTS="-Xmx${MEM_GB}G -Djava.io.tmpdir=./tmp"

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
  -h, --help                  Show this help message

Examples:
  $0 -s sample1
  $0 -s sample1 -t 16 -m 32
  $0 -s sample1 --r1 /path/to/sample1_R1.fastq.gz --r2 /path/to/sample1_R2.fastq.gz

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
            GATK_JAVA_OPTS="-Xmx${MEM_GB}G -Djava.io.tmpdir=./tmp"
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
DBSNP_VCF="${REF_DIR}/dbsnp_146.hg38.vcf.gz"
MILLS_INDELS_VCF="${REF_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
THOUSAND_GENOMES_VCF="${REF_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

# Create output directories
SAMPLE_RESULTS_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_RESULTS_DIR}"/{bam,qc,recal,gvcf,logs}
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

# Check required tools
log_info "Checking required tools..."
for cmd in fastqc bwa samtools gatk; do
    check_command "$cmd"
done

# Check input files
log_info "Checking input files..."
check_file "$FASTQ_R1"
check_file "$FASTQ_R2"
check_file "$REF_GENOME"
check_file "$EXOME_TARGETS"
check_file "$DBSNP_VCF"
check_file "$MILLS_INDELS_VCF"

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
    gatk --java-options "${GATK_JAVA_OPTS}" MarkDuplicates \
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

# --- STEP 4: Base Quality Score Recalibration (BQSR) ---
log_info "=== STEP 4: Base Quality Score Recalibration ==="
RECAL_TABLE="${SAMPLE_RESULTS_DIR}/recal/${SAMPLE_NAME}.recal_data.table"
ANALYSIS_READY_BAM="${SAMPLE_RESULTS_DIR}/bam/${SAMPLE_NAME}.analysis_ready.bam"

if [ ! -f "${ANALYSIS_READY_BAM}" ]; then
    # Step 4a: Build the recalibration model
    log_info "Building BQSR recalibration model..."
    gatk --java-options "${GATK_JAVA_OPTS}" BaseRecalibrator \
        -R "${REF_GENOME}" \
        -I "${MARKED_DUPS_BAM}" \
        --known-sites "${DBSNP_VCF}" \
        --known-sites "${MILLS_INDELS_VCF}" \
        --known-sites "${THOUSAND_GENOMES_VCF}" \
        -L "${EXOME_TARGETS}" \
        -O "${RECAL_TABLE}" \
        2> "${SAMPLE_RESULTS_DIR}/logs/base_recalibrator.log"

    # Step 4b: Apply the recalibration
    log_info "Applying BQSR recalibration..."
    gatk --java-options "${GATK_JAVA_OPTS}" ApplyBQSR \
        -R "${REF_GENOME}" \
        -I "${MARKED_DUPS_BAM}" \
        -bqsr "${RECAL_TABLE}" \
        -L "${EXOME_TARGETS}" \
        -O "${ANALYSIS_READY_BAM}" \
        2> "${SAMPLE_RESULTS_DIR}/logs/apply_bqsr.log"

    log_success "BQSR completed"
else
    log_info "Analysis-ready BAM already exists, skipping BQSR"
fi

# --- STEP 5: Variant Calling ---
log_info "=== STEP 5: Variant Calling ==="
GVCF_FILE="${SAMPLE_RESULTS_DIR}/gvcf/${SAMPLE_NAME}.g.vcf.gz"

if [ ! -f "${GVCF_FILE}" ]; then
    log_info "Calling variants with HaplotypeCaller in GVCF mode..."
    gatk --java-options "${GATK_JAVA_OPTS}" HaplotypeCaller \
        -R "${REF_GENOME}" \
        -I "${ANALYSIS_READY_BAM}" \
        -O "${GVCF_FILE}" \
        -L "${EXOME_TARGETS}" \
        -ERC GVCF \
        --native-pair-hmm-threads "${THREADS}" \
        2> "${SAMPLE_RESULTS_DIR}/logs/haplotype_caller.log"
    log_success "Variant calling completed"
else
    log_info "GVCF file already exists, skipping variant calling"
fi

# --- STEP 6: Collect Quality Metrics ---
log_info "=== STEP 6: Quality Metrics ==="
HS_METRICS="${SAMPLE_RESULTS_DIR}/qc/${SAMPLE_NAME}.hs_metrics.txt"

if [ ! -f "${HS_METRICS}" ]; then
    log_info "Collecting hybrid selection metrics..."
    gatk --java-options "${GATK_JAVA_OPTS}" CollectHsMetrics \
        -I "${ANALYSIS_READY_BAM}" \
        -O "${HS_METRICS}" \
        -R "${REF_GENOME}" \
        -BAIT_INTERVALS "${EXOME_TARGETS}" \
        -TARGET_INTERVALS "${EXOME_TARGETS}" \
        2> "${SAMPLE_RESULTS_DIR}/logs/collect_hs_metrics.log"
    log_success "Quality metrics collection completed"
else
    log_info "Quality metrics already collected, skipping"
fi

# --- STEP 7: Generate Summary Report ---
log_info "=== STEP 7: Generate Summary Report ==="
SUMMARY_REPORT="${SAMPLE_RESULTS_DIR}/${SAMPLE_NAME}_summary.txt"

cat > "${SUMMARY_REPORT}" << EOF
GATK WES Pipeline Summary Report
================================
Sample: ${SAMPLE_NAME}
Pipeline Version: ${SCRIPT_VERSION}
Analysis Date: $(date)
Reference Genome: $(basename "${REF_GENOME}")

Input Files:
- R1: ${FASTQ_R1}
- R2: ${FASTQ_R2}
- Exome Targets: ${EXOME_TARGETS}

Output Files:
- Analysis-ready BAM: ${ANALYSIS_READY_BAM}
- gVCF: ${GVCF_FILE}
- Quality Metrics: ${HS_METRICS}
- Duplicate Metrics: ${METRICS_FILE}

Processing Parameters:
- Threads: ${THREADS}
- Memory: ${MEM_GB}GB
- Java Options: ${GATK_JAVA_OPTS}

Read Group Information:
- ID: ${RG_ID}
- Sample: ${RG_SM}
- Library: ${RG_LB}
- Platform: ${RG_PL}
- Platform Unit: ${RG_PU}

Log Files:
- FastQC: ${SAMPLE_RESULTS_DIR}/logs/fastqc_*.log
- BWA-MEM: ${SAMPLE_RESULTS_DIR}/logs/bwa_mem.log
- Samtools: ${SAMPLE_RESULTS_DIR}/logs/samtools_sort.log
- Mark Duplicates: ${SAMPLE_RESULTS_DIR}/logs/mark_duplicates.log
- BQSR: ${SAMPLE_RESULTS_DIR}/logs/base_recalibrator.log, ${SAMPLE_RESULTS_DIR}/logs/apply_bqsr.log
- HaplotypeCaller: ${SAMPLE_RESULTS_DIR}/logs/haplotype_caller.log
- Quality Metrics: ${SAMPLE_RESULTS_DIR}/logs/collect_hs_metrics.log

Next Steps:
1. Review quality metrics and FastQC reports
2. For cohort analysis, collect all gVCF files and run joint genotyping
3. Use scripts/joint_genotyping.sh for cohort-level analysis

EOF

log_success "Summary report generated: ${SUMMARY_REPORT}"

# --- COMPLETION ---
log_success "=== PIPELINE COMPLETED SUCCESSFULLY ==="
log_info "Per-sample processing for ${SAMPLE_NAME} completed"
log_info "Total runtime: $SECONDS seconds"
log_info "Output gVCF: ${GVCF_FILE}"
log_info "Summary report: ${SUMMARY_REPORT}"
echo
log_info "To perform joint genotyping on multiple samples, use:"
log_info "  scripts/joint_genotyping.sh -i results/gatk/*/gvcf/*.g.vcf.gz"

log_success "GATK Pipeline finished successfully!"