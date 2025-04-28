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

```

---
config:
  theme: neo
  look: neo
  layout: dagre
---
flowchart TB
 subgraph Data_and_List_Generation["Data & List Generation"]
        PD["Pathogenic Databases"]
        H["Human sources: papers & pipelines (16SPIP, FAPROTAX, etc.)"]
        F["Fish sources: Wardeh et al."]
        A["Arthropod sources: Wardeh et al."]
        Cf["Crustacean filtering (ensembl_crustacea.py)"]
        PPL["Potential pathogenic species lists"]
        SE["Synonym expansion (get_sp_synonyms.py)"]
        EXT["Extended PPB species lists<br>(Human / Fish / Crustacea)"]
  end
 subgraph Stage1["Stage 1: Database Builder"]
        B1["1_run_database_builder.py"]
        INIT["Initial datasets"]
        UNF["Unfiltered FASTA per group"]
        RAW["Raw metadata"]
  end
 subgraph Stage2["Stage 2: DB Filters"]
        F2["2_run_db_filters.py"]
        FFA["Filtered FASTA"]
        INT["Intermediate files (optional)"]
  end
 subgraph Stage3["Stage 3: DB Curation"]
        C1["3_run_db_curation.py"]
        TAX["Tax_reviewed_(group).xlsx"]
        CFA["Unfiltered FASTA"]
        SUM["Species-level summary"]
  end
 subgraph s1["<span style=padding-left:>Stage 1b: BACTERIA Database </span>"]
        SILVA["SILVA SSU db (v138.2)"]
        BF["Filter to Bacteria only"]
        BFC@{ label: "Remove 'uncultured' & 'unidentified'" }
  end
    PD --> H & F & A
    H --> PPL
    A --> Cf
    Cf --> PPL
    SE --> EXT
    SILVA --> BF
    BF --> BFC
    EXT --> INIT
    INIT --> B1
    B1 --> UNF & RAW & s1
    UNF --> F2
    F2 --> FFA & INT
    FFA --> C1
    C1 --> TAX & CFA & SUM
    CFA --> OUT1["data/output/(group)/Curated_DB.fasta"]
    SUM --> OUT2["all_groups_general_results.csv"]
    BFC --> B1
    F --> PPL
    PPL --> SE
    TAX --> CFA
    BFC@{ shape: rect}
```
