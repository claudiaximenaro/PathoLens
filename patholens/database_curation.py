import os
from patholens.fasta_processing import filter_fasta_only_retained_taxonomy
from patholens.config import get_group_output_dir

def build_curated_database(species_to_exclude_file, group):
   

    # Get output directory for the group
    group_output_dir = get_group_output_dir(group)
    fasta_file=os.path.join(group_output_dir, f"{group}_Pathogen_DB_Unfiltered.fasta")

    filter_fasta_only_retained_taxonomy(species_to_exclude_file, fasta_file, group, group_output_dir)

    print(f"Curated Database for {group} successfully generated in {group_output_dir}.")
