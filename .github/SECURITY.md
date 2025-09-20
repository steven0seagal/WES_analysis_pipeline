# Security Policy

## Supported Versions

We currently support the following versions of the WES Analysis Pipeline with security updates:

| Version | Supported          |
| ------- | ------------------ |
| 1.x.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Reporting a Vulnerability

We take the security of the WES Analysis Pipeline seriously. If you discover a security vulnerability, please report it responsibly.

### How to Report

1. **Do not** create a public GitHub issue for security vulnerabilities
2. Send an email to [security@yourorganization.com] with:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact
   - Suggested fix (if any)

### What to Expect

- **Response time**: We aim to respond within 48 hours
- **Assessment**: We will assess the vulnerability and determine severity
- **Fix timeline**: Critical vulnerabilities will be patched within 7 days
- **Disclosure**: We follow responsible disclosure practices

### Security Considerations for Users

When using this pipeline:

1. **Input validation**: Always validate input FASTQ and reference files
2. **File permissions**: Ensure proper file permissions in shared environments
3. **Container security**: Keep Docker/Singularity images updated
4. **Network security**: Be cautious when downloading reference data
5. **Cluster security**: Follow your institution's HPC security guidelines

### Dependencies

This pipeline relies on external tools and databases. Users should:

- Keep bioinformatics tools updated
- Use official Docker images when available
- Verify checksums for downloaded reference data
- Review conda package sources

### Cluster Environments

When running on shared HPC systems:

- Use appropriate resource limits
- Secure temporary directories
- Clean up intermediate files
- Follow institutional data handling policies