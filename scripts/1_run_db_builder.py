import os
import argparse
from patholens.database_builder import build_database
from patholens.config import INPUT_DIR, SILVA_DIR

def main():
    parser = argparse.ArgumentParser(description="PathoLens - Pathogenic Bacteria Database.")
    
    parser.add_argument("--fasta", required=True, help="FASTA input file (SILVA).")
    parser.add_argument("--species", required=True, help="TXT file with pathogenic species.")
    parser.add_argument("--group", required=True, help="Taxonomic group (HUMAN, FISH, CRUSTACEAN).")

    args = parser.parse_args()

    fasta_file = os.path.join(SILVA_DIR, args.fasta)
    species_file = os.path.join(INPUT_DIR, args.group, args.species)
    
    build_database(fasta_file, species_file, args.group)

if __name__ == "__main__":
    main()
