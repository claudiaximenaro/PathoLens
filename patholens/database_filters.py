import os
import pandas as pd
from patholens.utils import save_to_txt, save_to_review_list
from patholens.filters import compare_taxonomy_genus_species, filter_gen_included
from patholens.config import get_group_output_dir

'''
Builds a curated taxonomy filter for a specific group by processing FASTA taxonomies,
applying genus-species comparison filters to identify discrepancies, genus-inclusion
filters to retain relevant taxonomies, and generating a review list for manual verification.
'''
def build_taxonomy_filters(fasta_file, group, save_intermediate=False):
   
    group_output_dir = get_group_output_dir(group)
    os.makedirs(group_output_dir, exist_ok=True)

    
    review_file_path = os.path.join(group_output_dir, f"Tax_to_manual-review_{group}.xlsx")

    # Step 1: Compare genus and species in the taxonomies
    discrepancies, complete_taxonomies = compare_taxonomy_genus_species(fasta_file,group)

    if save_intermediate:
        save_to_txt(os.path.join(group_output_dir, f"discrepancies_taxonomies_{group}.txt"), discrepancies)

    # Step 2: Filter taxonomies where genus is included
    gen_included, to_review = filter_gen_included(discrepancies, group, complete_taxonomies)

    if save_intermediate:
        save_to_txt(os.path.join(group_output_dir, f"gen_included_{group}.txt"), gen_included)
        save_to_txt(os.path.join(group_output_dir, f"to_review_{group}.txt"), to_review)

    save_to_review_list(to_review, review_file_path)

    print(f"Taxonomy filtering for {group} completed. Results saved in {group_output_dir}.")
