"""
pciSeq: Probabilistic Cell Typing for Spatial Transcriptomics

A Python package for analyzing spatial transcriptomics data with probabilistic cell typing.
This setup script handles package configuration and dependencies.

For development installation:
    pip install -e .

For building distribution:
    python -m build
"""

import os
from pathlib import Path
from setuptools import setup, find_packages


def get_static_files(root: str) -> list[str]:
    """
    Collect all static files needed for web interface.

    Args:
        root: Root directory to search for static files

    Returns:
        List of file paths relative to pciSeq package
    """
    static_extensions = {'.html', '.js', '.css', '.msi'}
    root_path = Path(root)

    static_files = []
    for path in root_path.rglob('*'):
        if path.suffix in static_extensions:
            relative_path = str(path.relative_to('pciSeq'))
            static_files.append(relative_path)

    return static_files


# Core dependencies required for basic functionality
INSTALL_REQUIRES = [
    'altair',           # Data visualization
    'dask',             # Parallel computing
    'diplib',           # Image processing
    'fastremap',        # Fast array operations
    'flask',            # Web framework
    'natsort',          # Natural sorting
    'numexpr',          # Fast numerical expressions
    'numpy_groupies',   # Group operations
    'pandas',           # Data manipulation
    'pyvips',           # Image processing
    'redis',            # Caching
    'scikit-image',     # Image analysis
    'scikit-learn',     # Machine learning
    'scipy',            # Scientific computing
    'streamlit',        # Web apps
    'tomlkit',          # TOML parsing
    'tqdm',             # Progress bars
    'colorlog',         # Colored logging
]

# Optional dependencies for interactive use
EXTRAS_REQUIRE = {
    'interactive': [
        'matplotlib>=2.2.0',
        'jupyter'
    ],
}


# Get version from _version.py
def get_version() -> str:
    """Extract version from _version.py."""
    version_file = Path('pciSeq/_version.py')
    if not version_file.exists():
        raise RuntimeError('Version file not found')

    version = None
    with open(version_file) as f:
        for line in f:
            if line.startswith('__version__'):
                version = line.split('=')[1].strip().strip("'")
                break

    if version is None:
        raise RuntimeError('Version not found in _version.py')

    return version


# Read long description from README
def get_long_description() -> str:
    """Read long description from README.md."""
    with open("README.md", encoding="utf-8") as f:
        return f.read()


setup(
    name="pciSeq",
    version=get_version(),
    license="BSD",
    author="Dimitris Nicoloutsopoulos",
    author_email="dimitris.nicoloutsopoulos@gmail.com",
    description="Probabilistic cell typing for spatial transcriptomics",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/acycliq/pciSeq",
    packages=find_packages(),
    python_requires=">=3.8",  # Specify minimum Python version
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    include_package_data=True,
    package_data={
        'pciSeq': get_static_files(os.path.join('pciSeq', 'static'))
    },
    classifiers=[
        "Development Status :: 4 - Beta",  # Add development status
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Image Processing",
    ],
    keywords=[
        "spatial-transcriptomics",
        "cell-typing",
        "bioinformatics",
        "image-analysis",
        "single-cell"
    ],
    project_urls={
        "Bug Reports": "https://github.com/acycliq/pciSeq/issues",
        "Source": "https://github.com/acycliq/pciSeq",
        "Documentation": "https://github.com/acycliq/pciSeq#readme"
    }
)
