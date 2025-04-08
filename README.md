# PathoLens

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

## Usage

### Running the database builder

```bash
python3 scripts/run_db_builder.py -h
```

#### Command-line arguments:

```
usage: run_db_builder.py [-h] --fasta FASTA --species SPECIES
                         --group GROUP [--unmatched UNMATCHED]

PathoLens - Pathogenic Bacteria Database.

options:
  -h, --help            show this help message and exit
  --fasta FASTA         FASTA input file (SILVA).
  --species SPECIES     TXT file with pathogenic species.
  --group GROUP         Taxonomic group (HUMAN, FISH, CRUSTACEAN).
  --unmatched UNMATCHED Optional: CSV file with additional
                        information for missing species table.
```

### Running the filters

```bash
python3 scripts/run_filters_db.py
```

### Running the curated database process

```bash
python3 scripts/run_curated_db.py
```

## Input Data Structure

PathoLens expects the following directory structure for input data:

```
data/
├── input/
│   ├── SILVA/
│   │   ├── SILVA_138.2_SSURef_tax_silva.fasta
│   ├── HUMAN/
│   │   ├── Human_patho_species.txt
│   │   ├── Human_not_silva.csv
│   ├── FISH/
│   ├── CRUSTACEAN/
```

## Important:

* You **do not need to provide full paths** for the input files, only the file names. The script assumes that the files are placed within the **structured input directories** imposed by PathoLens.

## Output Data Structure

Results are saved in the following directory structure:

```
data/
├── output/
│   ├── HUMAN/
│   │   ├── HUMAN_Pathogen_DB.fasta
│   ├── FISH/
│   ├── CRUSTACEAN/
│   ├── all_groups_general_results.csv
```

## Configuration

PathoLens uses a `config.py` file to manage input and output directories.

## License

This project is licensed under the MIT License.

## Contact

For questions or contributions, please contact claudia.restrepo-ortiz@ird.fr or open an issue on GitHub.
