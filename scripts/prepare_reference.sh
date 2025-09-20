#!/bin/bash

# prepare_reference.sh
# This script downloads the GRCh38 reference genome and related resources,
# then prepares all necessary indices for BWA, Samtools, and GATK.
# It also downloads SnpEff databases for variant annotation.

set -euo pipefail  # Exit on error, unset variable, or pipe failure

# --- Configuration ---
REF_DIR="reference"
GENOME_FASTA="GRCh38.p13.genome.fa"
DBSNP_VCF="dbsnp_146.hg38.vcf.gz"
MILLS_INDELS_VCF="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
THOUSAND_GENOMES_VCF="1000G_phase1.snps.high_confidence.hg38.vcf.gz"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

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
    if ! command -v "$1" &> /dev/null; then
        log_error "Command '$1' not found. Please ensure it's installed and in PATH."
        exit 1
    fi
}

download_with_retry() {
    local url="$1"
    local output="$2"
    local max_attempts=3
    local attempt=1

    while [ $attempt -le $max_attempts ]; do
        log_info "Download attempt $attempt/$max_attempts for $(basename "$output")"
        if wget -c -O "$output" "$url"; then
            log_success "Successfully downloaded $(basename "$output")"
            return 0
        else
            log_warning "Download attempt $attempt failed"
            ((attempt++))
            sleep 5
        fi
    done

    log_error "Failed to download $output after $max_attempts attempts"
    return 1
}

# --- Main Script ---

log_info "Starting reference data preparation for WES Analysis Pipeline"

# Check required commands
log_info "Checking required tools..."
for cmd in wget bwa samtools gatk java; do
    check_command "$cmd"
done
log_success "All required tools found"

# Create reference directory
log_info "Creating reference directory: $REF_DIR"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# --- Download Reference Genome ---
log_info "Downloading GRCh38 reference genome..."
if [ ! -f "$GENOME_FASTA" ]; then
    download_with_retry \
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz" \
        "${GENOME_FASTA}.gz"

    log_info "Decompressing reference genome..."
    gunzip "${GENOME_FASTA}.gz"
    log_success "Reference genome ready: $GENOME_FASTA"
else
    log_info "Reference genome already exists, skipping download"
fi

# --- Download Known Variant Sites ---
log_info "Downloading known variant sites for GATK BQSR..."

# dbSNP
if [ ! -f "$DBSNP_VCF" ]; then
    download_with_retry \
        "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \
        "$DBSNP_VCF"
    download_with_retry \
        "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi" \
        "${DBSNP_VCF}.tbi"
else
    log_info "dbSNP VCF already exists, skipping download"
fi

# Mills and 1000G Gold Standard Indels
if [ ! -f "$MILLS_INDELS_VCF" ]; then
    download_with_retry \
        "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
        "$MILLS_INDELS_VCF"
    download_with_retry \
        "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
        "${MILLS_INDELS_VCF}.tbi"
else
    log_info "Mills indels VCF already exists, skipping download"
fi

# 1000 Genomes high confidence SNPs
if [ ! -f "$THOUSAND_GENOMES_VCF" ]; then
    download_with_retry \
        "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
        "$THOUSAND_GENOMES_VCF"
    download_with_retry \
        "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" \
        "${THOUSAND_GENOMES_VCF}.tbi"
else
    log_info "1000 Genomes VCF already exists, skipping download"
fi

# --- Create Indices ---
log_info "Creating genome indices..."

# BWA index
if [ ! -f "${GENOME_FASTA}.bwt" ]; then
    log_info "Creating BWA index..."
    bwa index "$GENOME_FASTA"
    log_success "BWA index created"
else
    log_info "BWA index already exists, skipping"
fi

# Samtools FASTA index
if [ ! -f "${GENOME_FASTA}.fai" ]; then
    log_info "Creating Samtools FASTA index..."
    samtools faidx "$GENOME_FASTA"
    log_success "Samtools FASTA index created"
else
    log_info "Samtools FASTA index already exists, skipping"
fi

# GATK Sequence Dictionary
DICT_FILE="${GENOME_FASTA%.fa}.dict"
if [ ! -f "$DICT_FILE" ]; then
    log_info "Creating GATK Sequence Dictionary..."
    gatk CreateSequenceDictionary -R "$GENOME_FASTA" -O "$DICT_FILE"
    log_success "GATK Sequence Dictionary created"
else
    log_info "GATK Sequence Dictionary already exists, skipping"
fi

# --- Download SnpEff Database ---
log_info "Setting up SnpEff annotation database..."
if command -v snpEff.jar &> /dev/null; then
    SNPEFF_CMD="snpEff.jar"
elif command -v snpEff &> /dev/null; then
    SNPEFF_CMD="snpEff"
else
    log_warning "SnpEff not found in PATH. You may need to download the database manually."
    log_info "Use: java -jar /path/to/snpEff.jar download -v GRCh38.86"
    SNPEFF_CMD=""
fi

if [ -n "$SNPEFF_CMD" ]; then
    log_info "Downloading SnpEff database for GRCh38..."
    if [ "$SNPEFF_CMD" = "snpEff.jar" ]; then
        java -jar "$(which snpEff.jar)" download -v GRCh38.86 || log_warning "SnpEff database download failed"
    else
        snpEff download -v GRCh38.86 || log_warning "SnpEff database download failed"
    fi
fi

# --- Create example exome targets file ---
log_info "Creating example exome targets file..."
cat > exome_targets.bed << 'EOF'
# This is an example exome targets file
# Users should replace this with their actual capture kit coordinates
# Format: chromosome start end [name] [score] [strand]
chr1	65419	65433
chr1	65520	65573
chr1	69037	70008
chr1	367659	368597
chr1	621059	622053
# Add more regions as needed for your specific capture kit
EOF

log_warning "IMPORTANT: Replace 'exome_targets.bed' with your actual capture kit coordinates!"

# --- Summary ---
cd ..  # Return to main directory

log_success "Reference data preparation completed!"
echo
log_info "Summary of downloaded/created files:"
echo "  ✓ Reference genome: reference/$GENOME_FASTA"
echo "  ✓ BWA index files: reference/${GENOME_FASTA}.*"
echo "  ✓ Samtools index: reference/${GENOME_FASTA}.fai"
echo "  ✓ GATK dictionary: reference/${DICT_FILE}"
echo "  ✓ dbSNP variants: reference/$DBSNP_VCF"
echo "  ✓ Mills indels: reference/$MILLS_INDELS_VCF"
echo "  ✓ 1000G SNPs: reference/$THOUSAND_GENOMES_VCF"
echo "  ✓ Example targets: reference/exome_targets.bed"
echo
log_warning "Remember to:"
echo "  1. Replace the example exome_targets.bed with your capture kit coordinates"
echo "  2. Ensure you have sufficient disk space (~30-50 GB for reference files)"
echo "  3. The SnpEff database may require additional manual setup"
echo
log_success "You can now run the WES analysis pipelines!"