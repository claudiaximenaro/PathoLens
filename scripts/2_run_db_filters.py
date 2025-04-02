import os
import argparse
from patholens.database_filters import build_taxonomy_filters
from patholens.config import get_group_output_dir

def main():
    parser = argparse.ArgumentParser(description="Run taxonomy filters on a FASTA file.")
    
    parser.add_argument("--group", required=True, help="Taxonomic group (HUMAN, FISH, CRUSTACEAN).")
    parser.add_argument("--save_intermediate", action='store_true', help="Save intermediate results.")

    args = parser.parse_args()

    group_output_dir = get_group_output_dir(args.group)
    os.makedirs(group_output_dir, exist_ok=True)

    fasta_file = os.path.join(group_output_dir, f"{args.group}_Pathogen_DB_Unfiltered.fasta")
    
    build_taxonomy_filters(fasta_file, args.group, args.save_intermediate)

if __name__ == "__main__":
    main()
