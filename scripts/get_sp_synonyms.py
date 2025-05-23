# This script fetches accepted scientific names and synonyms for a list of species from NCBI Entrez. 
# It processes species data and outputs the results in an Excel file, including synonyms and taxonomic information.



import argparse
import pandas as pd
from Bio import Entrez
import time
import re

def clean_species_name(full_name):
    """Extract only the binomial name (Genus + species), discard additional details, and remove brackets."""
    full_name = re.sub(r"\[([A-Za-z]+)\]", r"\1", full_name).strip()
    match = re.match(r"([A-Z][a-z]+) ([a-z]+)", full_name)
    return f"{match.group(1)} {match.group(2)}" if match else full_name

def get_species_info(species_name, email, max_retries=3):
    Entrez.email = email
    retries = 0
    """Get the TaxID and accepted scientific names for a species."""
    while retries < max_retries:
        try:
            clean_name = clean_species_name(species_name)  
            handle = Entrez.esearch(db="taxonomy", term=clean_name)
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                print(f"No tax ID found for {clean_name}. Trying with AccNumber...")
                return None, clean_name, []

            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            tax_record = Entrez.read(handle)[0]

            accepted_name = tax_record.get("ScientificName", clean_name)
            other_names = tax_record.get("OtherNames", {})

            synonyms = other_names.get("Synonym", [])
            includes = other_names.get("Includes", [])
            basionym = other_names.get("Basionym", [])

            all_synonyms = synonyms + basionym + includes
            return tax_id, accepted_name, all_synonyms
        except Exception as e:
            retries += 1
            print(f"Error fetching data for {species_name} (attempt {retries}/{max_retries}): {e}")
            time.sleep(2 ** retries) 
            return None, species_name, []

def get_valid_accession_number(species_name, max_results=5):
    """Search for an AccNumber in nuccore and obtain a valid organism."""
    try:
        search_handle = Entrez.esearch(db="nuccore", term=species_name, retmax=max_results, retmode="xml")
        search_record = Entrez.read(search_handle)

        if not search_record["IdList"]:
            print(f"No accession number found for {species_name}")
            return None, None

        for acc_id in search_record["IdList"]:
            fetch_handle = Entrez.efetch(db="nuccore", id=acc_id, rettype="gb", retmode="xml")
            fetch_record = Entrez.read(fetch_handle)

            for entry in fetch_record:
                if "GBSeq_organism" in entry:
                    organism = entry["GBSeq_organism"]
                    if "unidentified" not in organism.lower():
                        clean_organism = clean_species_name(organism)
                        return acc_id, clean_organism
                    else:
                        print(f"Ignoring {organism} (ID: {acc_id}) - unidentified.")

        return None, None
    except Exception as e:
        print(f"Error fetching AccNumber for {species_name}: {e}")
        return None, None

def generate_species_extended_list (synonyms_excel, group):
    df = pd.read_excel(synonyms_excel, header=0)
    df = df.loc[:, df.columns.notna()]
    sp_name = df.values.flatten()
    sp_name = [str(valor).strip() for valor in sp_name if pd.notna(valor) and str(valor).strip()]

    uniques_sp = list(dict.fromkeys(sp_name))

    with open(f'data/input/{group}/{group}_sp_pathogens_list-EXT.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(uniques_sp))  


def process_species_file(input_file, email, output_file):

    """Process the species list, fetching the accepted name and synonyms."""
    with open(input_file, 'r', encoding='utf-8') as file:
        species_list = [line.strip() for line in file if line.strip()]
    
    results = []
    max_synonyms = 0

    for species in species_list:
        time.sleep(0.5)  # Avoid overloading the NCBI API
        tax_id, accepted_name, synonyms = get_species_info(species, email)

        if tax_id is None:  # If no TaxID was found, try with AccNumber
            acc_number, organism = get_valid_accession_number(species)
            if organism:
                print(f"Obtained {organism} via AccNumber")
                tax_id, accepted_name, synonyms = get_species_info(organism, email)

        eid_value = species if species != accepted_name else ""
        row = [accepted_name, eid_value] + synonyms
        max_synonyms = max(max_synonyms, len(synonyms))
        results.append(row)

    headers = ["Accepted Scientific Name", "EID_sp"] + [f"Synonym {i+1}" for i in range(max_synonyms)]
    
    for row in results:
        row.extend([""] * (max_synonyms - (len(row) - 2)))
    
    df = pd.DataFrame(results, columns=headers)
    
    
    df = df.sort_values(by=["EID_sp"], ascending=False).drop_duplicates(subset=["Accepted Scientific Name"], keep='first')
    df = df.sort_values(by=["Accepted Scientific Name"])

    df.to_excel(output_file, index=False)
    print(f"File saved: {output_file}")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch species synonyms from NCBI and save to Excel.")
    parser.add_argument("input_file", help="Path to the input text file containing species names.")
    parser.add_argument("email", help="Email address required by NCBI Entrez API.")
    parser.add_argument("group", help="Taxonomic group for the species.")
    
    args = parser.parse_args()
    output_file=f'data/input/{args.group}/{args.group}_Pathogen_TaxSyn_List.xlsx'
    process_species_file(args.input_file, args.email, output_file)
    print("..........")
    generate_species_extended_list (output_file, args.group)
    print("Process completed.")
