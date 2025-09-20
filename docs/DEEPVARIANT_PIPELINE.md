# DeepVariant Pipeline Documentation

This document provides detailed information about the DeepVariant pipeline implementation for Whole Exome Sequencing analysis.

## Overview

The DeepVariant pipeline uses Google's DeepVariant algorithm for variant calling. DeepVariant employs a deep neural network to call genetic variants from aligned reads, providing high accuracy for both SNPs and indels. This pipeline is particularly effective for WES data and offers an alternative to traditional GATK-based approaches.

## Pipeline Steps

### 1. Quality Control (FastQC)

**Purpose**: Assess the quality of raw sequencing reads

**Tools**: FastQC v0.12.1

**Input**: Raw FASTQ files

**Output**: Quality control reports

**Same as GATK pipeline** - refer to GATK documentation for details.

### 2. Read Alignment (BWA-MEM)

**Purpose**: Align reads to the reference genome

**Tools**: BWA v0.7.17

**Input**: FASTQ files + reference genome

**Output**: Sorted BAM file

**Same parameters as GATK pipeline** - refer to GATK documentation.

### 3. Mark Duplicates (Picard)

**Purpose**: Identify and mark duplicate reads

**Tools**: Picard v3.1.1

**Input**: Sorted BAM file

**Output**: Duplicate-marked BAM

**Same as GATK pipeline** - no BQSR step needed for DeepVariant.

### 4. Variant Calling (DeepVariant)

**Purpose**: Call variants using deep learning

**Tools**: DeepVariant v1.6.1 (via Docker/Singularity)

**Input**: Analysis-ready BAM

**Output**: gVCF file

**Key Features**:
- Deep neural network for variant calling
- High accuracy for SNPs and indels
- Containerized execution
- Optimized for WES data

**Parameters**:
- Model: WES-specific model
- Number of shards: 8-16
- Regions: Exome targets

### 5. Annotation (SnpEff)

**Purpose**: Add functional annotations

**Tools**: SnpEff v5.2

**Input**: VCF file

**Output**: Annotated VCF

**Same as GATK pipeline**.

## DeepVariant Specific Features

### Container Execution

DeepVariant runs in a container for reproducibility:

```bash
# Docker execution
docker run \
    -v /path/to/data:/data \
    gcr.io/deepvariant-docker/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=/data/reference.fa \
    --reads=/data/sample.bam \
    --regions=/data/targets.bed \
    --output_vcf=/data/output.vcf.gz \
    --output_gvcf=/data/output.g.vcf.gz \
    --num_shards=8
```

### Model Selection

- **WES Model**: Optimized for whole exome sequencing
- **WGS Model**: For whole genome sequencing
- **PacBio Model**: For long-read data

### Performance Characteristics

- **Accuracy**: Higher sensitivity for indels
- **Speed**: Faster than GATK for large datasets
- **Memory**: Lower memory requirements
- **Scalability**: Better parallel processing

## Configuration Parameters

### DeepVariant Parameters

```yaml
# config/config.yaml
deepvariant:
  version: "1.6.1"
  model_type: "WES"
  num_shards: 16
  container_runtime: "docker"
```

### Container Settings

```yaml
advanced:
  use_deepvariant_container: true
  container_runtime: "docker"
```

## Quality Control Metrics

### DeepVariant Specific Metrics

- **Variant Quality Score**: DeepVariant's internal quality metric
- **Allele Frequency**: Population frequency estimates
- **Read Support**: Number of reads supporting variant
- **Mapping Quality**: Average mapping quality of supporting reads

### Expected Performance

- **SNP Sensitivity**: >99%
- **SNP Precision**: >99.9%
- **Indel Sensitivity**: >95%
- **Indel Precision**: >99%

## Comparison with GATK

### Strengths of DeepVariant

1. **Higher Indel Accuracy**: Better detection of insertions/deletions
2. **Consistent Performance**: Less sensitive to coverage variations
3. **Faster Processing**: Reduced computational requirements
4. **Simplified Workflow**: No BQSR step required

### When to Use DeepVariant

- High indel burden expected
- Variable coverage across targets
- Limited computational resources
- Research requiring high indel sensitivity

### When to Use GATK

- Well-established protocols required
- Extensive customization needed
- Integration with existing GATK workflows
- Regulatory compliance requiring GATK

## Troubleshooting

### Container Issues

#### Docker Permission Denied

```bash
# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker
```

#### Container Pull Failures

```bash
# Pull manually
docker pull gcr.io/deepvariant-docker/deepvariant:1.6.1

# Check network connectivity
docker run hello-world
```

### Model Loading Issues

**Symptoms**: Model fails to load

**Solutions**:
- Verify model type matches data
- Check container version compatibility
- Ensure sufficient memory available

### Output Validation

**Symptoms**: Empty or invalid VCF

**Solutions**:
- Check input BAM validity
- Verify regions file format
- Review container logs

## Performance Optimization

### Sharding Strategy

```bash
# Optimal shard count
num_shards = min(available_cores, num_targets / 1000)
```

### Memory Management

- **Per Shard**: ~2-4 GB RAM
- **Total Memory**: num_shards × 4GB + overhead

### Parallel Processing

```bash
# Run multiple samples in parallel
parallel -j 4 'bash deepvariant_pipeline.sh -s {}' ::: sample1 sample2 sample3 sample4
```

## Integration Considerations

### Joint Calling with GLnexus

DeepVariant gVCFs can be joint-called using GLnexus:

```bash
# Install GLnexus
conda install -c bioconda glnexus

# Joint calling
glnexus \
    --config DeepVariantWES \
    sample1.g.vcf.gz \
    sample2.g.vcf.gz \
    > cohort.vcf
```

### Comparison with GATK Joint Genotyping

- **GLnexus**: Faster, simpler
- **GATK GenotypeGVCFs**: More configurable, better for large cohorts

## Best Practices

### Data Preparation

1. **BAM Sorting**: Ensure BAM is coordinate-sorted
2. **Index Validation**: Verify BAM index is current
3. **Region Filtering**: Use exome targets for efficiency

### Quality Assurance

1. **Positive Controls**: Include reference samples
2. **Cross-Validation**: Compare results with GATK
3. **Reproducibility**: Use fixed random seeds if available

### Monitoring

1. **Resource Usage**: Track CPU and memory per shard
2. **Runtime**: Monitor processing time per sample
3. **Output Size**: Check VCF file sizes for reasonableness

## Future Developments

### DeepVariant v2.0

- Improved model architecture
- Better indel calling
- Enhanced performance
- New model types

### Integration Improvements

- Native Snakemake support
- Better cluster integration
- Automated model selection

This pipeline provides a modern, accurate approach to WES variant calling using deep learning, offering advantages in speed and indel detection compared to traditional methods.