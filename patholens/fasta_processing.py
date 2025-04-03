from Bio import SeqIO
import os
from collections import Counter
from patholens.utils import load_species, count_and_extract_taxonomies, extract_species_to_remove,save_to_txt

def analyze_columns(file_path, categories=None):
    """    
    Parameters:
        file_path (str): Path to the FASTA file.
        categories (dict, optional): A dictionary where keys are category names and values are
                                     the prefix strings to identify the taxonomy.
                                     Defaults to {"Bacteria": "Bacteria", "Eukaryota": "Eukaryota"}.
    
    Returns:
        tuple: A tuple containing:
            - taxonomy_counts (Counter): Counts of sequences keyed by the number of taxonomic columns.
            - category_examples (dict): For each category, a dict mapping column count to an example header.
            - category_counts (dict): Total counts of sequences for each category.
    """
    if categories is None:
        categories = {"Bacteria": "Bacteria", "Eukaryota": "Eukaryota"}
        
    taxonomy_counts = Counter()
    category_examples = {cat: {} for cat in categories}
    category_counts = {cat: 0 for cat in categories}
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line.strip()
                # Extract taxonomy part after the first space (if available)
                taxonomy = header.split(' ', 1)[1] if ' ' in header else ""
                columns = taxonomy.split(';')
                num_columns = len(columns)
                taxonomy_counts[num_columns] += 1
                
                # Check each category for a match based on the taxonomy prefix
                for cat, prefix in categories.items():
                    if taxonomy.startswith(prefix):
                        category_counts[cat] += 1
                        if num_columns not in category_examples[cat]:
                            category_examples[cat][num_columns] = header
                        
    return taxonomy_counts, category_examples, category_counts



def filter_fasta(input_fasta, output_fasta, filter_func):
    """
    Filters a FASTA file using a provided filtering function.
    
    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file.
        filter_func (callable): A function that takes a SeqRecord and returns True
                                if the record should be included, False otherwise.
    """
    # Parse all records from the input FASTA file
    records = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Apply the filtering function to each record
    filtered_records = [record for record in records if filter_func(record)]
    
    # Write the filtered records to the output FASTA file
    with open(output_fasta, "w") as out_handle:
        SeqIO.write(filtered_records, out_handle, "fasta")

# Define a filtering function for species-based filtering
def species_filter(record):
    """
    Returns True if the record's species (extracted as the last part of the description
    after splitting by ';') does not contain 'uncultured' or 'unidentified'.
    """
    species = record.description.split(';')[-1].strip().lower()
    return not any(term in species for term in ('uncultured', 'unidentified'))

# Define a factory to create a filtering function based on taxonomy prefix
def taxonomy_filter(taxonomy_prefix):
    """
    Returns a filtering function that checks if the taxonomy (the part of the header
    after the first space) starts with the given taxonomy_prefix.
    """
    def filter_func(record):
        # Extract taxonomy part after the first space (if available)
        taxonomy = record.description.split(' ', 1)[1] if ' ' in record.description else ""
        return taxonomy.startswith(taxonomy_prefix)
    return filter_func


def extract_fasta_group(fasta_file, species_file, output_fasta, group):
    """Extract sequences from a FASTA file that match species in the species list."""
    species_set = load_species(species_file, return_type="set")
    
    with open(output_fasta, 'w', encoding='utf-8') as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            description = record.description.split(";")   
            species_name = description[-1]
            
            # Check if any species name from the list is in the description
            if species_name in species_set:
                
                out_fasta.write(f">{record.description}\n{str(record.seq)}\n")
    
    count_and_extract_taxonomies(output_fasta, group,"Initial")

def filter_fasta_only_retained_taxonomy(exclude_taxa_file, fasta_input, group, group_output,save_intermediate):
    """
    Filters a FASTA file by removing sequences with exactly matching taxonomies.

    Parameters:
    - taxonomies (set): Set of full taxonomies to remove.
    - fasta_input (str): Path to the input FASTA file to filter.
    - fasta_output (str): Path to the output FASTA file after filtering.

    Returns:
    - tuple: (number of sequences removed, number of unique taxonomies removed)
    """
    exclude_taxa=extract_species_to_remove(exclude_taxa_file)
    # Define output file paths inside the group output folder
    fasta_output = os.path.join(group_output, f"{group}_Pathogen_DB_curated.fasta")


    with open(fasta_input, 'r') as f:
        fasta_lines = f.readlines()

    filtered_fasta = []
    exclude = False
    removed_sequences = 0
    removed_taxonomies = set()

    for line in fasta_lines:
        if line.startswith('>'):  
            parts = line.strip().split(' ', 1)  
            if len(parts) > 1:
                taxonomy = parts[1]  
                exclude = taxonomy in exclude_taxa
                if exclude:
                    removed_sequences += 1
                    removed_taxonomies.add(taxonomy)
        if not exclude:
            filtered_fasta.append(line)

    with open(fasta_output, 'w') as f:
        f.writelines(filtered_fasta)
    
    if save_intermediate:
        uniques_tax,_=count_and_extract_taxonomies(fasta_output, group, "Curated")
        save_to_txt(os.path.join(group_output, f"final_unique_taxonomies_{group}.txt"), uniques_tax)
    else:
        count_and_extract_taxonomies(fasta_output, group, "Curated")

    return removed_sequences, len(removed_taxonomies)
