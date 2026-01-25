# PathoLens 🔬🦠

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

> Extract and curate pathogenic bacteria databases from SILVA for HUMAN, FISH, and CRUSTACEAN hosts.

## Overview

PathoLens is a Python package designed to extract and curate a database of pathogenic bacteria from the SILVA database. It filters sequences based on taxonomic criteria and generates a structured output for further analysis.

## Features

* Extracts only bacterial sequences from the SILVA database.
* Filters out undefined and uncultured species.
* Processes pathogenic species lists for different taxonomic groups (HUMAN, FISH, CRUSTACEAN).
* Generates intermediate and final curated FASTA files and metadata reports.

## Installation

Clone the repository and install the package along with its dependencies:

```bash
 git clone https://github.com/claudiaximenaro/PathoLens.git
 cd PathoLens
 pip install -e .
```

## Requirements

- Python >= 3.8
- biopython >= 1.79
- pandas >= 1.3.0
- openpyxl >= 3.0.0
- requests >= 2.25.0

All dependencies are automatically installed when running `pip install -e .`

## Usage

## Input Data Structure

PathoLens expects the following directory structure for input data, make sure you have at least the SILVA database and the list of species in their respective folders:

```
.
├── data
│   ├── input
│   │   ├── CRUSTACEAN
│   │   │   ├── Crustacean_sp_pathogens_list.txt
│   │   ├── FISH
│   │   ├── HUMAN
│   │   └── SILVA
│   │       └── SILVA_138.2_SSURef_tax_silva.fasta
```

## Important:

* You **do not need to provide full paths** for the input files, only the file names. The script assumes that the files are placed within the **structured input directories** imposed by PathoLens.

## Quick Start Example

```bash
# Step 1: Build the database for HUMAN pathogens
python3 scripts/1_run_db_builder.py \
    --fasta SILVA_138.2_SSURef_tax_silva.fasta \
    --species Human_sp_pathogens_list.txt \
    --group HUMAN

# Step 2: Apply taxonomy filters
python3 scripts/2_run_db_filters.py --group HUMAN --save_intermediate

# Step 3: Curate the final database
python3 scripts/3_run_db_curation.py \
    --sp_remove Tax_to_manual-review_HUMAN.xlsx \
    --group HUMAN
```

After installation with `pip install -e .`, you can also use the CLI commands:

```bash
patholens-build --fasta SILVA_138.2_SSURef_tax_silva.fasta --species Human_sp_pathogens_list.txt --group HUMAN
patholens-filter --group HUMAN --save_intermediate
patholens-curate --sp_remove Tax_to_manual-review_HUMAN.xlsx --group HUMAN
```

## Pipeline Workflow

```
┌──────────────────┐     ┌──────────────────┐     ┌──────────────────┐
│  SILVA Database  │     │  Species List    │     │  Exclusion File  │
│  (FASTA)         │     │  (TXT)           │     │  (Excel)         │
└────────┬─────────┘     └────────┬─────────┘     └────────┬─────────┘
         │                        │                        │
         ▼                        ▼                        │
    ┌─────────────────────────────────────┐                │
    │      1. Database Builder            │                │
    │   Filters bacteria + pathogens      │                │
    └─────────────────┬───────────────────┘                │
                      ▼                                    │
    ┌─────────────────────────────────────┐                │
    │      2. Taxonomy Filters            │                │
    │   Removes discrepancies             │                │
    └─────────────────┬───────────────────┘                │
                      ▼                                    │
    ┌─────────────────────────────────────┐                │
    │      3. Database Curation           │◄───────────────┘
    │   Final curated database            │
    └─────────────────┬───────────────────┘
                      ▼
    ┌─────────────────────────────────────┐
    │     OUTPUT: Curated FASTA DB        │
    │  + Species reports (CSV/Excel)      │
    └─────────────────────────────────────┘
```

  
## Command-line arguments:

### Running the database builder

```bash
python3 scripts/1_run_db_builder.py -h
```
```
usage: 1_run_db_builder.py [-h] --fasta FASTA --species SPECIES
                           --group GROUP

PathoLens - Pathogenic Bacteria Database.

options:
  -h, --help         show this help message and exit
  --fasta FASTA      FASTA input file (SILVA).
  --species SPECIES  TXT file with pathogenic species.
  --group GROUP      Taxonomic group (HUMAN, FISH, CRUSTACEAN).
```

### Running the filters

```bash
python3 scripts/2_run_db_filters.py -h
```
```
usage: 2_run_db_filters.py [-h] --group GROUP [--save_intermediate]

Run taxonomy filters on a FASTA file.

options:
  -h, --help           show this help message and exit
  --group GROUP        Taxonomic group (HUMAN, FISH, CRUSTACEAN).
  --save_intermediate  Save intermediate results. (optional)
```
### Running the curated database process

```bash
python3 scripts/3_run_db_curation.py -h
```
```
usage: 3_run_db_curation.py [-h] --sp_remove SP_REMOVE --group GROUP

Curated Pathogen Database Generator.

options:
  -h, --help            show this help message and exit
  --sp_remove SP_REMOVE Excel file with species to exclude.                   
  --group GROUP         Taxonomic group (HUMAN, FISH, CRUSTACEAN).
```

## Output Data Structure

Results are saved in the following directory structure:

```
.
├── data
│   ├── input
│   ├── output
│   │   ├── CRUSTACEAN
│   │   ├── FISH
│   │   ├── HUMAN
│   │       └── HUMAN_Pathogen_DB.fasta
│   └── all_groups_general_results.csv

```

## Configuration

PathoLens uses a `config.py` file to manage input and output directories.

## License

This project is licensed under the GPL-3.0 license.

## Contact

For questions or contributions, please contact claudia.restrepo-ortiz@ird.fr or open an issue on GitHub.
