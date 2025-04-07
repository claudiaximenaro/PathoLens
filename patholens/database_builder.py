import os
from patholens.fasta_processing import filter_fasta, species_filter, taxonomy_filter, extract_fasta_group
from patholens.utils import extract_and_filter_species, find_missing_species, process_files, count_and_extract_taxonomies
from patholens.config import SILVA_DIR, get_group_output_dir

'''Builds a curated FASTA database for a specific group of pathogens 
by filtering and processing sequences from the SILVA database, applying
taxonomic filters to retain only bacterial sequences, species-level filters
to remove unidentified or uncultured species, and group-specific filters 
to extract pathogen sequences relevant to the target group'''


def build_database(fasta_file, species_file, group):
    group_output_dir = get_group_output_dir(group)


    bacteria_fasta = os.path.join(SILVA_DIR, "Bacteria.fasta")
    filtered_fasta = os.path.join(SILVA_DIR, "Bacteria_filtered.fasta")
    group_path_fasta = os.path.join(group_output_dir, f"{group}_Pathogen_DB_Unfiltered.fasta")
    
    

    
    if not os.path.exists(bacteria_fasta) or not os.path.exists(filtered_fasta):
        print("Bacteria FASTA files not found. Generating them...")
        filter_fasta(fasta_file, bacteria_fasta, taxonomy_filter("Bacteria"))  # Filter only bacteria
        filter_fasta(bacteria_fasta, filtered_fasta, species_filter)  # Remove unidentified species
        print(f"Bacteria FASTA files generated and stored in {SILVA_DIR}")
    else:
        print("Skipping construction of Bacteria-only FASTA. Using existing files.")

    count_and_extract_taxonomies(fasta_file, "SILVA", "Initial")
    count_and_extract_taxonomies(bacteria_fasta, "SILVA", "Bacteria")
    count_and_extract_taxonomies(filtered_fasta, "Bacteria", "uncultured-unidentified")
        
   
    extract_fasta_group(filtered_fasta, species_file, group_path_fasta, group)

    extract_and_filter_species(group_path_fasta, species_file, group)

    find_missing_species(group_path_fasta, species_file, group)

    print(f"Database for {group} successfully generated in {group_output_dir}.")
