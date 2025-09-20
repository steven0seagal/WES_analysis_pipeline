#!/usr/bin/env python3

"""
Setup script for WES Analysis Pipeline
"""

from setuptools import setup, find_packages
import os

# Read README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
def read_requirements(filename):
    with open(filename, "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]

# Package metadata
setup(
    name="wes-pipeline",
    version="1.0.0",
    author="WES Pipeline Team",
    author_email="wes@example.com",
    description="Whole Exome Sequencing Analysis Pipeline with GATK and DeepVariant",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/wes-pipeline",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    keywords="bioinformatics genomics wes sequencing gatk deepvariant",
    python_requires=">=3.9",
    install_requires=read_requirements("requirements.txt"),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "pre-commit",
        ],
        "docs": [
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
    entry_points={
        "console_scripts": [
            "wes-pipeline=wes_pipeline.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "wes_pipeline": [
            "config/*.yaml",
            "scripts/*.sh",
            "docs/*",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/yourusername/wes-pipeline/issues",
        "Source": "https://github.com/yourusername/wes-pipeline",
        "Documentation": "https://wes-pipeline.readthedocs.io/",
    },
)