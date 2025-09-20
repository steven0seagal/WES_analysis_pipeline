# Contributing to WES Analysis Pipeline

Thank you for your interest in contributing to the WES Analysis Pipeline! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Code Standards](#code-standards)
- [Testing](#testing)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)
- [Reporting Issues](#reporting-issues)

## Code of Conduct

This project follows a code of conduct to ensure a welcoming environment for all contributors. By participating, you agree to:

- Be respectful and inclusive
- Focus on constructive feedback
- Accept responsibility for mistakes
- Show empathy towards other contributors
- Help create a positive community

## Getting Started

### Prerequisites

- Git
- Conda/Mamba
- Docker (for DeepVariant)
- Python 3.9+
- Shell scripting knowledge

### Development Setup

1. Fork the repository on GitHub
2. Clone your fork:
```bash
git clone https://github.com/yourusername/wes-pipeline.git
cd wes-pipeline
```

3. Set up the development environment:
```bash
conda env create -f environment.yml
conda activate wes-analysis
```

4. Install development dependencies:
```bash
pip install -r requirements-dev.txt
```

5. Set up pre-commit hooks:
```bash
pre-commit install
```

## Development Workflow

### Branching Strategy

- `main`: Production-ready code
- `develop`: Integration branch for features
- `feature/*`: Feature branches
- `bugfix/*`: Bug fix branches
- `hotfix/*`: Critical fixes for production

### Creating a Feature Branch

```bash
git checkout develop
git pull origin develop
git checkout -b feature/your-feature-name
```

### Making Changes

1. Write tests for new functionality
2. Implement the feature
3. Run tests and linting
4. Update documentation
5. Commit changes

## Code Standards

### Shell Scripts

- Use `#!/bin/bash`
- Include `set -euo pipefail`
- Use descriptive variable names
- Include error handling
- Add comments for complex logic
- Follow POSIX shell standards where possible

### Python Code

- Follow PEP 8 style guide
- Use type hints
- Write docstrings for functions and classes
- Use meaningful variable names
- Keep functions small and focused

### YAML Configuration

- Use consistent indentation (2 spaces)
- Add comments for complex configurations
- Use descriptive keys
- Validate syntax

### Commit Messages

Follow conventional commit format:

```
type(scope): description

[optional body]

[optional footer]
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `style`: Code style changes
- `refactor`: Code refactoring
- `test`: Testing
- `chore`: Maintenance

Examples:
```
feat(pipeline): add support for custom reference genomes
fix(alignment): resolve BWA memory issue with large genomes
docs(readme): update installation instructions
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_gatk_pipeline.py

# Run with coverage
pytest --cov=src --cov-report=html
```

### Test Structure

- Unit tests for individual functions
- Integration tests for pipeline components
- End-to-end tests for complete workflows
- Performance tests for benchmarking

### Test Data

- Use small, representative test datasets
- Include positive and negative controls
- Test edge cases and error conditions
- Validate outputs against expected results

## Documentation

### Documentation Types

- **README.md**: Project overview and quick start
- **INSTALL.md**: Detailed installation guide
- **USAGE.md**: Usage instructions and examples
- **API Documentation**: Function and class documentation
- **Technical Docs**: Detailed technical specifications

### Documentation Standards

- Use Markdown format
- Include code examples
- Provide troubleshooting guides
- Keep documentation up-to-date
- Use clear, concise language

## Submitting Changes

### Pull Request Process

1. Ensure all tests pass
2. Update documentation
3. Run linting and formatting
4. Create a pull request
5. Request review from maintainers

### Pull Request Template

Please use the following template:

```markdown
## Description
Brief description of the changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Refactoring
- [ ] Testing

## Testing
- [ ] Unit tests added/updated
- [ ] Integration tests added/updated
- [ ] Manual testing performed

## Checklist
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] Tests pass
- [ ] No breaking changes
```

### Code Review

- All pull requests require review
- Address review comments promptly
- Maintain respectful communication
- Focus on code quality and functionality

## Reporting Issues

### Bug Reports

When reporting bugs, please include:

- **Description**: Clear description of the issue
- **Steps to Reproduce**: Step-by-step instructions
- **Expected Behavior**: What should happen
- **Actual Behavior**: What actually happens
- **Environment**: OS, versions, hardware
- **Logs**: Error messages and relevant logs

### Feature Requests

For feature requests, please include:

- **Description**: Clear description of the feature
- **Use Case**: Why this feature is needed
- **Implementation Ideas**: Suggestions for implementation
- **Alternatives**: Alternative approaches considered

### Issue Labels

- `bug`: Something isn't working
- `enhancement`: New feature or improvement
- `documentation`: Documentation issues
- `question`: Questions or discussions
- `help wanted`: Good for newcomers
- `good first issue`: Simple issues for beginners

## Recognition

Contributors will be recognized in:
- GitHub repository contributors
- Changelog entries
- Release notes
- Project documentation

## Getting Help

- **Documentation**: Check the docs/ directory
- **Issues**: Search existing issues on GitHub
- **Discussions**: Use GitHub Discussions for questions
- **Community**: Join our community forum

Thank you for contributing to the WES Analysis Pipeline!