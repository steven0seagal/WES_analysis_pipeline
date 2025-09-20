#!/bin/bash

# lint_all.sh
# Run all linting and quality checks on the codebase

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo "Running comprehensive code quality checks..."
echo "Project root: $PROJECT_ROOT"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

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

cd "$PROJECT_ROOT"

# Check if required tools are installed
check_tools() {
    log_info "Checking required tools..."

    local missing_tools=()

    if ! command -v shellcheck &> /dev/null; then
        missing_tools+=("shellcheck")
    fi

    if ! command -v black &> /dev/null; then
        missing_tools+=("black")
    fi

    if ! command -v flake8 &> /dev/null; then
        missing_tools+=("flake8")
    fi

    if ! command -v yamllint &> /dev/null; then
        missing_tools+=("yamllint")
    fi

    if [ ${#missing_tools[@]} -ne 0 ]; then
        log_warning "Missing tools: ${missing_tools[*]}"
        log_info "Install with: pip install black flake8 yamllint"
        log_info "Install shellcheck with: apt-get install shellcheck (Ubuntu/Debian)"
        return 1
    fi

    log_success "All required tools are installed"
    return 0
}

# Lint shell scripts
lint_shell() {
    log_info "Linting shell scripts..."

    local shell_files
    shell_files=$(find scripts/ -name "*.sh" -type f)

    if [ -z "$shell_files" ]; then
        log_warning "No shell scripts found"
        return 0
    fi

    local failed=0
    while IFS= read -r file; do
        if ! shellcheck "$file"; then
            log_error "Shellcheck failed for: $file"
            failed=1
        fi
    done <<< "$shell_files"

    if [ $failed -eq 0 ]; then
        log_success "All shell scripts passed linting"
    fi

    return $failed
}

# Lint Python files
lint_python() {
    log_info "Linting Python files..."

    local python_files
    python_files=$(find . -name "*.py" -type f -not -path "./.*")

    if [ -z "$python_files" ]; then
        log_warning "No Python files found"
        return 0
    fi

    local failed=0

    # Check formatting with black
    log_info "Checking Python formatting with black..."
    if ! black --check --diff .; then
        log_error "Black formatting check failed"
        log_info "Run 'black .' to fix formatting"
        failed=1
    else
        log_success "Python formatting is correct"
    fi

    # Lint with flake8
    log_info "Linting Python with flake8..."
    if ! flake8 .; then
        log_error "Flake8 linting failed"
        failed=1
    else
        log_success "Python linting passed"
    fi

    return $failed
}

# Lint YAML files
lint_yaml() {
    log_info "Linting YAML files..."

    local yaml_files
    yaml_files=$(find . -name "*.yml" -o -name "*.yaml" | grep -v "\.github")

    if [ -z "$yaml_files" ]; then
        log_warning "No YAML files found"
        return 0
    fi

    local failed=0
    while IFS= read -r file; do
        if ! yamllint "$file"; then
            log_error "YAML linting failed for: $file"
            failed=1
        fi
    done <<< "$yaml_files"

    if [ $failed -eq 0 ]; then
        log_success "All YAML files passed linting"
    fi

    return $failed
}

# Check file permissions
check_permissions() {
    log_info "Checking file permissions..."

    local failed=0

    # Check that shell scripts are executable
    local non_executable
    non_executable=$(find scripts/ -name "*.sh" -type f ! -executable)

    if [ -n "$non_executable" ]; then
        log_error "Non-executable shell scripts found:"
        echo "$non_executable"
        log_info "Run: chmod +x scripts/*.sh"
        failed=1
    else
        log_success "All shell scripts are executable"
    fi

    return $failed
}

# Check for large files
check_large_files() {
    log_info "Checking for large files..."

    local large_files
    large_files=$(find . -type f -size +50M -not -path "./.*" -not -path "./reference/*" -not -path "./results/*")

    if [ -n "$large_files" ]; then
        log_warning "Large files found (should be in .gitignore):"
        echo "$large_files"
        log_info "Consider adding large files to .gitignore"
    else
        log_success "No large files found in tracked directories"
    fi
}

# Main execution
main() {
    local overall_failed=0

    if ! check_tools; then
        log_error "Required tools missing. Install dependencies and try again."
        exit 1
    fi

    if ! lint_shell; then
        overall_failed=1
    fi

    if ! lint_python; then
        overall_failed=1
    fi

    if ! lint_yaml; then
        overall_failed=1
    fi

    if ! check_permissions; then
        overall_failed=1
    fi

    check_large_files

    if [ $overall_failed -eq 0 ]; then
        log_success "All quality checks passed!"
        exit 0
    else
        log_error "Some quality checks failed. Please fix the issues above."
        exit 1
    fi
}

# Run main function
main "$@"