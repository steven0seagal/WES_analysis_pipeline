# Changelog

All notable changes to the WES Analysis Pipeline will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial release of WES Analysis Pipeline
- GATK Best Practices pipeline implementation
- DeepVariant pipeline implementation
- Snakemake workflow orchestration
- Docker/Singularity container support
- Comprehensive documentation
- CI/CD with GitHub Actions
- Pre-commit hooks for code quality
- Automated testing framework

### Changed
- N/A (initial release)

### Deprecated
- N/A (initial release)

### Removed
- N/A (initial release)

### Fixed
- N/A (initial release)

### Security
- N/A (initial release)

## [1.0.0] - 2024-01-XX

### Added
- Complete GATK Best Practices workflow
  - Quality control with FastQC
  - BWA-MEM alignment and sorting
  - Picard MarkDuplicates
  - Base Quality Score Recalibration (BQSR)
  - HaplotypeCaller in GVCF mode
  - Hybrid selection metrics collection
- DeepVariant pipeline implementation
  - Container-based execution
  - Optimized for WES data
  - GLnexus joint calling support
- Snakemake workflows for both pipelines
- Comprehensive configuration system
- Reference data preparation scripts
- Quality control and validation tools
- MultiQC integration for report generation
- SnpEff annotation support

### Documentation
- Complete README with installation and usage
- Detailed installation guide
- Usage instructions and examples
- Pipeline-specific technical documentation
- Troubleshooting guides
- API documentation

### CI/CD
- GitHub Actions workflows for CI
- Automated testing with small datasets
- Code quality checks (linting, formatting)
- Pre-commit hooks configuration
- Docker image building

### Development
- Conda environment specification
- Development dependencies
- Testing framework setup
- Code quality tools integration

---

## Types of Changes

- `Added` for new features
- `Changed` for changes in existing functionality
- `Deprecated` for soon-to-be removed features
- `Removed` for now removed features
- `Fixed` for any bug fixes
- `Security` in case of vulnerabilities

## Versioning

This project uses [Semantic Versioning](https://semver.org/):

- **MAJOR** version for incompatible API changes
- **MINOR** version for backwards-compatible functionality additions
- **PATCH** version for backwards-compatible bug fixes

## Release Process

1. Update version numbers in relevant files
2. Update CHANGELOG.md with release notes
3. Create git tag for the release
4. Push to main branch
5. Create GitHub release
6. Update documentation if needed

## Contributing to Changelog

- Keep entries brief but descriptive
- Group similar changes together
- Use present tense for changes
- Reference issue numbers when applicable
- Update the Unreleased section during development