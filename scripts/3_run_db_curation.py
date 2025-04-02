import os
import argparse
from patholens.database_curation import build_curated_database

def main():
    parser = argparse.ArgumentParser(description="Curated Pathogen Database Generator.")
    
    parser.add_argument("--sp_remove", required=True, help="Excel file with species to exclude.")
    parser.add_argument("--group", required=True, help="Taxonomic group (HUMAN, FISH, CRUSTACEAN).")
    parser.add_argument("--save_intermediate", action='store_true', help="Save intermediate results.")

    args = parser.parse_args()

    build_curated_database(args.sp_remove, args.group, args.save_intermediate)

if __name__ == "__main__":
    main()

