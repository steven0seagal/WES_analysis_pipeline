# Contributing to WES Analysis Pipeline

We welcome contributions! This document provides guidelines for contributing to the WES Analysis Pipeline project.

## Getting Started

### Prerequisites

- Git
- Conda/Miniconda
- Docker (for DeepVariant testing)
- Basic knowledge of bioinformatics workflows

### Development Setup

1. **Fork and clone the repository:**
   ```bash
   git clone https://github.com/your-username/wes-analysis-pipeline.git
   cd wes-analysis-pipeline
   ```

2. **Create the development environment:**
   ```bash
   conda env create -f environment.yml
   conda activate wes-analysis
   ```

3. **Install pre-commit hooks:**
   ```bash
   pre-commit install
   ```

4. **Test your setup:**
   ```bash
   ./scripts/test_environment_locally.sh
   ```

## How to Contribute

### Reporting Issues

- Use GitHub Issues to report bugs or request features
- Search existing issues before creating new ones
- Provide detailed information including:
  - Operating system and version
  - Conda/Python version
  - Steps to reproduce the issue
  - Expected vs actual behavior
  - Relevant log files

### Contributing Code

#### Development Workflow

1. **Create a feature branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes:**
   - Follow coding standards
   - Add tests for new functionality
   - Update documentation as needed

3. **Test locally:**
   ```bash
   ./scripts/test_environment_locally.sh
   ./scripts/validate_installation.sh
   pre-commit run --all-files
   ```

4. **Commit your changes:**
   ```bash
   git add .
   git commit -m "feat: add new feature description"
   ```

5. **Push and create a pull request:**
   ```bash
   git push origin feature/your-feature-name
   ```

## Coding Standards

### Shell Scripts
- Use `#!/bin/bash` shebang
- Enable strict mode: `set -euo pipefail`
- Add comments for complex logic
- Follow shellcheck recommendations

### Python Code
- Follow PEP 8 style guide
- Use Black for code formatting
- Maximum line length: 88 characters

### Snakemake Workflows
- Use descriptive rule names
- Include proper input validation
- Follow Snakemake best practices

## License

By contributing, you agree that your contributions will be licensed under the MIT License.