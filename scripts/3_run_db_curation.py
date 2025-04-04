import os
import argparse
from patholens.database_curation import build_curated_database
from patholens.database_analysis import extract_final_species
from patholens.config import get_group_output_dir


def main():
    parser = argparse.ArgumentParser(description="Curated Pathogen Database Generator.")
    
    parser.add_argument("--sp_remove", required=True, help="Excel file with species to exclude.")
    parser.add_argument("--group", required=True, help="Taxonomic group (HUMAN, FISH, CRUSTACEAN).")

    args = parser.parse_args()

    build_curated_database(args.sp_remove, args.group)
    
    group_output_dir = get_group_output_dir(args.group)
    os.makedirs(group_output_dir, exist_ok=True)
    file_outpath = os.path.join(group_output_dir, f"Species_match_{args.group}.xlsx")

    extract_final_species(file_outpath,args.group,group_output_dir)

if __name__ == "__main__":
    main()

