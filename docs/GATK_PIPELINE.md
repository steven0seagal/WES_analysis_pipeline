# GATK Pipeline Documentation

This document provides detailed information about the GATK Best Practices pipeline implementation for Whole Exome Sequencing analysis.

## Overview

The GATK pipeline follows the [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) workflow for germline short variant discovery. It processes raw FASTQ files through alignment, preprocessing, and variant calling to produce high-quality gVCF files suitable for joint genotyping.

## Pipeline Steps

### 1. Quality Control (FastQC)

**Purpose**: Assess the quality of raw sequencing reads

**Tools**: FastQC v0.12.1

**Input**: Raw FASTQ files (`sample_R1.fastq.gz`, `sample_R2.fastq.gz`)

**Output**:
- `qc/sample_R1_fastqc.html` - Quality report for R1
- `qc/sample_R2_fastqc.html` - Quality report for R2
- `qc/sample_R1_fastqc.zip` - Raw data for R1
- `qc/sample_R2_fastqc.zip` - Raw data for R2

**Key Metrics to Check**:
- Per base sequence quality
- Per sequence quality scores
- Per base sequence content
- Adapter content
- Overrepresented sequences

### 2. Read Alignment (BWA-MEM)

**Purpose**: Align reads to the reference genome

**Tools**: BWA v0.7.17

**Input**: FASTQ files + reference genome

**Output**: `bam/sample.sorted.bam` - Coordinate-sorted BAM file

**Parameters**:
- Algorithm: BWA-MEM (maximal exact matches)
- Read group information included
- Duplicate reads marked for removal

**Command**:
```bash
bwa mem -t 16 -R "@RG\\tID:group1.sample\\tPL:ILLUMINA\\tPU:unit1.sample\\tLB:lib1.sample\\tSM:sample\\tCN:SequencingCenter\\tDS:WES" reference/GRCh38.p13.genome.fa sample_R1.fastq.gz sample_R2.fastq.gz | samtools sort -@ 16 -o bam/sample.sorted.bam -
```

### 3. Mark Duplicates (Picard)

**Purpose**: Identify and mark duplicate reads

**Tools**: Picard v3.1.1

**Input**: Sorted BAM file

**Output**:
- `bam/sample.marked_dups.bam` - BAM with duplicates marked
- `qc/sample.dup_metrics.txt` - Duplication metrics

**Key Metrics**:
- **Duplication Rate**: Percentage of reads marked as duplicates
- **Library Size**: Estimated library size based on duplication
- **Optical Duplicates**: Duplicates due to optical clustering

**Acceptable Values**:
- Duplication rate: <20% for exome sequencing
- Optical duplicates: <5% of total duplicates

### 4. Base Quality Score Recalibration (BQSR)

**Purpose**: Correct systematic errors in base quality scores

**Tools**: GATK v4.5.0.0

**Input**: Duplicate-marked BAM + known variant sites

**Output**:
- `recal/sample.recal_data.table` - Recalibration table
- `bam/sample.analysis_ready.bam` - Recalibrated BAM

**Known Sites Used**:
- dbSNP v146 (common variants)
- Mills indels (known indels)
- 1000 Genomes indels (additional indels)

**Steps**:
1. **BaseRecalibrator**: Build recalibration model
2. **ApplyBQSR**: Apply recalibration to reads

### 5. Variant Calling (HaplotypeCaller)

**Purpose**: Call germline variants in GVCF format

**Tools**: GATK HaplotypeCaller

**Input**: Analysis-ready BAM

**Output**: `gvcf/sample.g.vcf.gz` - Genomic VCF

**Parameters**:
- Mode: GVCF (genomic VCF)
- Emit reference confidence: GVCF
- Target regions: Exome targets BED file
- Native pair HMM threads: 8

**Key Features**:
- Joint calling compatible
- Reference confidence scores
- Handles indels and SNPs
- Parallel processing support

### 6. Quality Metrics Collection

**Purpose**: Collect comprehensive quality metrics

**Tools**: GATK CollectHsMetrics

**Input**: Analysis-ready BAM

**Output**: `qc/sample.hs_metrics.txt` - Hybrid selection metrics

**Key Metrics**:
- **Mean Target Coverage**: Average coverage across targets
- **Coverage Uniformity**: Percentage of targets at various coverage levels
- **Target Coverage**: Percentage of targets covered at 20x, 50x, etc.
- **Fold Enrichment**: On-target vs off-target coverage ratio

### 7. Annotation (SnpEff)

**Purpose**: Add functional annotations to variants

**Tools**: SnpEff v5.2

**Input**: VCF file from joint genotyping

**Output**: Annotated VCF with functional predictions

**Databases**:
- GRCh38.86 (Ensembl)
- Includes: gene names, functional class, impact, etc.

## Configuration Parameters

### Computational Resources

```yaml
# config/config.yaml
params:
  threads: 16
  mem_gb: 32
  java_opts: "-Xmx32G -Djava.io.tmpdir=./tmp"
```

### GATK-Specific Parameters

```yaml
gatk:
  base_recalibrator_threads: 8
  haplotype_caller_threads: 8
  joint_genotyping_threads: 16
```

### Quality Thresholds

- **Mapping Quality**: MQ ≥ 20
- **Base Quality**: Q ≥ 20
- **Depth**: DP ≥ 10
- **Genotype Quality**: GQ ≥ 20

## Performance Considerations

### Memory Requirements

- **BaseRecalibrator**: 4-8 GB per thread
- **HaplotypeCaller**: 2-4 GB per thread
- **Joint Genotyping**: 8-16 GB per thread

### Runtime Estimates

- **Small Exome** (30-50Mb): 2-4 hours
- **Large Exome** (50-100Mb): 4-8 hours
- **Whole Genome**: 24-48 hours

### Optimization Strategies

1. **Parallel Processing**: Use multiple threads for I/O intensive steps
2. **Memory Management**: Adjust Java heap size based on available RAM
3. **Disk I/O**: Use SSD storage for reference data and temporary files
4. **Targeted Analysis**: Limit processing to exome targets

## Quality Control Checks

### Pre-Analysis

- [ ] FASTQ files exist and are not corrupted
- [ ] Reference genome and indices are present
- [ ] Known variant files are indexed
- [ ] Exome targets BED file is valid
- [ ] Sufficient disk space available

### Post-Analysis

- [ ] BAM files are coordinate-sorted and indexed
- [ ] Duplicate rate is within acceptable range
- [ ] Coverage metrics meet requirements
- [ ] gVCF files are valid and contain variants
- [ ] No errors in pipeline logs

## Troubleshooting

### Common Issues

#### High Duplication Rate

**Symptoms**: Duplication rate >30%

**Possible Causes**:
- Over-amplification during library preparation
- Optical duplicates from sequencing
- PCR artifacts

**Solutions**:
- Review library preparation protocol
- Adjust Picard MarkDuplicates parameters
- Consider molecular barcoding

#### Low Coverage Uniformity

**Symptoms**: Poor coverage distribution across targets

**Possible Causes**:
- GC bias in library preparation
- Capture efficiency issues
- Sequencing artifacts

**Solutions**:
- Review capture kit performance
- Adjust alignment parameters
- Consider bait redesign

#### Variant Calling Artifacts

**Symptoms**: Unexpected variant patterns

**Possible Causes**:
- Reference genome issues
- BQSR overfitting
- Mapping artifacts

**Solutions**:
- Verify reference genome version
- Check BQSR known sites
- Review mapping quality filters

## Best Practices

### Data Management

1. **Backup Raw Data**: Always keep original FASTQ files
2. **Version Control**: Track pipeline versions and parameters
3. **Documentation**: Record all analysis parameters and decisions
4. **Validation**: Cross-validate results with alternative methods

### Quality Assurance

1. **Positive Controls**: Include known variant samples
2. **Replicate Analysis**: Run samples in duplicate
3. **Batch Effects**: Monitor for technical artifacts
4. **Reference Materials**: Use NIST or similar reference materials

### Performance Monitoring

1. **Resource Usage**: Track CPU, memory, and disk usage
2. **Runtime Metrics**: Monitor step completion times
3. **Error Rates**: Track pipeline failure rates
4. **Quality Metrics**: Establish baseline quality thresholds

## Integration with Downstream Analysis

### Joint Genotyping

After individual sample processing:

```bash
# Combine gVCF files
gatk GenomicsDBImport \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    --genomicsdb-workspace-path workspace \
    -L exome_targets.bed

# Joint genotyping
gatk GenotypeGVCFs \
    -R reference.fa \
    -V gendb://workspace \
    -O cohort.vcf.gz
```

### Variant Filtering

Apply hard filters or use VQSR:

```bash
# Hard filtering
gatk VariantFiltration \
    -V cohort.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8.0" \
    -O cohort.filtered.vcf.gz
```

### Functional Annotation

```bash
# SnpEff annotation
snpEff -v GRCh38.86 cohort.filtered.vcf.gz > cohort.annotated.vcf

# Extract fields
gatk VariantsToTable \
    -V cohort.annotated.vcf \
    -F CHROM -F POS -F REF -F ALT -F ANN \
    -O variant_table.txt
```

This pipeline provides a robust foundation for WES analysis, following GATK best practices while maintaining flexibility for customization and optimization.