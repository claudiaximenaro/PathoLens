
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from patholens.config import OUTPUT_DIR, GENERAL_RESULTS_FILE

def load_species(file_path, return_type="set", verbose=False):
    """
    Loads species names from a text file.
    
    Parameters:
        file_path (str): Path to the species file.
        return_type (str): The type of collection to return: "set" (default) or "list".
        verbose (bool): If True, prints the number of species loaded.
        
    Returns:
        A set or a list of species names, based on the return_type parameter.
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        species = [line.strip() for line in f if line.strip()]
    
    if return_type.lower() == "set":
        species = set(species)
    
    if verbose:
        print(f"Loaded {len(species)} species.")
    
    return species


def write_to_csv_general_results(total_sequences, unique_taxonomies_count, group, filter_step):
    """ save statistics. """
    if not group or not filter_step:
        print("Group and filter step are required to save results.")
        return

    if os.path.exists(GENERAL_RESULTS_FILE):
        df = pd.read_csv(GENERAL_RESULTS_FILE, sep='\t')

        if ((df["Group"] == group) & (df["Filter"] == filter_step)).any():
            print(f"Results for {group} - {filter_step} already exist in CSV. Skipping...")
            return

    new_data = {
        "Group": group,
        "Filter": filter_step,
        "N° of Sequences": total_sequences,
        "N° of Unique Taxonomies": unique_taxonomies_count
    }

    new_row = pd.DataFrame([new_data])

    # Guardar en el CSV
    if os.path.exists(GENERAL_RESULTS_FILE):
        new_row.to_csv(GENERAL_RESULTS_FILE, mode='a', sep='\t', index=False, header=False)
    else:
        new_row.to_csv(GENERAL_RESULTS_FILE, sep='\t', index=False)

    print(f"Saved results for {group} - {filter_step} in CSV.")


def count_and_extract_taxonomies(fasta_file, group=None, filter_step=None):
    total_sequences = 0
    unique_descriptions = set()
    sequences = defaultdict(int)
    complete_taxonomies = []

    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                total_sequences += 1
                description = ' '.join(line.strip().split()[1:])  
                unique_descriptions.add(description)
                complete_taxonomies.append(description)
                sequences[description] += 1

    print(f"Total sequences in {group}-->{filter_step}: {total_sequences}")
    print(f"Unique descriptions in {group}-->{filter_step}: {len(unique_descriptions)}")

    
    if group and filter_step:
        write_to_csv_general_results(total_sequences, len(unique_descriptions), group, filter_step)

    return unique_descriptions, complete_taxonomies

def save_to_txt(filename, data):
    with open(filename, 'w') as f:
        for item in data:
            f.write(f"{item}\n")
    print(f"Saved {len(data)} items to {filename}")


def save_to_review_list(to_review, output_file_path):
    df_to_review = pd.DataFrame(to_review, columns=["Complete Taxonomy"])
    df_to_review["Retained"] = ""  
    df_to_review.to_excel(output_file_path, index=False)
    print(f"Saved to-review list to {output_file_path}")


def check_synonyms(actual_names, unique_species_list):
    """Check if any of the actual names exist in the unique species list."""
    return any(name in unique_species_list for name in actual_names)



def process_files(missing_species_file, unmatched_file, unique_species_file, group):
    """
    Processes species and unmatched data files, filtering and matching species, then saves the results.

    Parameters:
        missing_species_file (str): Path to the file containing missing species.
        unmatched_file (str): Path to the unmatched species file.
        unique_species_file (str): Path to the file containing unique species.
        group (str): Group name used for file naming.

    Output:
        - "<group>_missing_sp_info.csv": Processed file with species filtering and matching.
    """

    output_file = os.path.join(OUTPUT_DIR, group, f"{group}_missing_sp_info.csv")

    missing_species_list = load_species(missing_species_file, return_type="set")
    unmatched_df = pd.read_csv(unmatched_file, sep=';', dtype=str)
    unique_species_list = load_species(unique_species_file, return_type="set")

    # Filter unmatched species file to keep only rows where 'sp_name' is in missing_species_list
    matched_df = unmatched_df[unmatched_df['sp_name'].isin(missing_species_list)].copy()

    # Check for synonyms in unique_species_list
    matched_df['found_by_synonym'] = matched_df.apply(
        lambda row: check_synonyms(
            [row['actual_sp_name_NCBI'], row['actual_sp_name_SILVA']]
            if pd.notna(row['actual_sp_name_NCBI']) and pd.notna(row['actual_sp_name_SILVA'])
            else [],
            unique_species_list
        ), axis=1
    )

    non_empty_columns = [col for col in matched_df.columns if col != 'sp_name']
    final_df = matched_df.dropna(how='all', subset=non_empty_columns)

    last_cols = [col for col in ['reason', 'comments'] if col in final_df.columns]
    cols = [col for col in final_df.columns if col not in last_cols] + last_cols
    final_df = final_df[cols]

    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed file saved to: {output_file}")



def find_missing_species(fasta_file, species_file, group):
    """
    EFinds species from the species file that are not present in the FASTA file.

    Parameters:
        fasta_file (str): Reference FASTA file.
        species_file (str): File containing the list of species.
        group (str): Prefix for the output files.

    Output:
        - "<group>_missing_species.txt": List of species missing from the FASTA file.
    """
    output_file = os.path.join(OUTPUT_DIR,group, f"{group}_missing_species.txt")

    species_set = load_species(species_file, return_type="set")
    found_species = set()


    for record in SeqIO.parse(fasta_file, "fasta"):
        description = record.description  
        for species in species_set:
            if species in description:
                found_species.add(species)

    missing_species = species_set - found_species
    print(f"Total species listed: {len(species_set)}")
    print(f"Total species found in {group}: {len(found_species)}")
    print(f"Total missing species: {len(missing_species)}")

    with open(output_file, 'w', encoding='utf-8') as out_f:
        for species in sorted(missing_species):
            out_f.write(species + "\n")

    print(f"List of missing species stored in: {output_file}")



def extract_and_filter_species(fasta_file, exclude_species_file, group):
    """
    Extracts unique species from a FASTA file, saves them to a unique species file, 
    and filters out species present in the exclusion file.
    
    Parameters:
        fasta_file (str): Path to the FASTA file.
        exclude_species_file (str): Path to a text file with species names to exclude (one per line).
        group (str): Group name to prefix the output files.
    
    Output:
        - "<group>_unique_species.txt" (contains all unique species)
        - "<group>_filtered_species.txt" (contains unique species after filtering)
    """
    
    unique_species_file = os.path.join(OUTPUT_DIR,group, f"{group}_unique_species.txt")
    filtered_species_file = os.path.join(OUTPUT_DIR,group, f"{group}_filtered_species.txt")

    species_set = set()
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                species = line.strip().split(';')[-1].strip()
                species_set.add(species)
    
    with open(unique_species_file, 'w') as out_file:
        for species in sorted(species_set):
            out_file.write(species + '\n')

    with open(exclude_species_file, 'r') as file:
        exclude_species = {line.strip() for line in file if line.strip()}
    
    filtered_species = species_set - exclude_species
    
    with open(filtered_species_file, 'w') as out_file:
        for species in sorted(filtered_species):
            out_file.write(species + '\n')

    print(f"Unique species saved to: {unique_species_file}")
    print(f"Filtered species saved to: {filtered_species_file}")

def extract_species_to_remove(removal_sp_file_path):
    """
    Extracts species marked for removal from the specified Excel file.

    Parameters:
    - file_path (str): Path to the Excel file containing species information.
    
    Returns:
    - list: A list of complete taxonomies to remove.
    """
    df = pd.read_excel(removal_sp_file_path)

    filtered_df = df[df['Retained'] == 'No']

    species_to_remove = filtered_df['Complete Taxonomy'].tolist()  
    return species_to_remove


def load_species_from_txt(txt_file):
    species_list = []
    with open(txt_file, 'r', encoding='utf-8') as file:
        for line in file:
            parts = line.strip().split(';')  
            if parts:
                species_list.append(parts[-1]) 
    return species_list


def load_excel_data(excel_file, sheet_name=0):
    return pd.read_excel(excel_file, sheet_name=sheet_name, dtype=str)

def find_species_in_excel(species_list, df):
    """Check if species are present in the Accepted Name or any Synonym column."""
    results = []
    
    
    df = df.fillna('')
   
    df['All Names'] = df.apply(lambda row: {row['Accepted Scientific Name'], row['EID_sp']} | set(row.iloc[2:].values), axis=1)
    
    for species in species_list:
        found = False
        for _, row in df.iterrows():
            if species in row['All Names']:
                results.append((species, row['Accepted Scientific Name']))
                found = True
                break
        if not found:
            results.append((species, 'Not Found'))
    
    return results