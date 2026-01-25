from setuptools import setup, find_packages

# Read the long description from README.md
try:
    with open("README.md", encoding="utf-8") as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "Extract and curate pathogenic bacteria databases from SILVA"

setup(
    name="patholens",
    version="0.1.0",
    author="Claudia Restrepo-Ortiz",
    author_email="claudia.restrepo-ortiz@ird.fr",
    description="Extract and curate pathogenic bacteria databases from SILVA",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/claudiaximenaro/PathoLens",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "openpyxl>=3.0.0",
        "requests>=2.25.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "patholens-build=scripts.1_run_db_builder:main",
            "patholens-filter=scripts.2_run_db_filters:main",
            "patholens-curate=scripts.3_run_db_curation:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)