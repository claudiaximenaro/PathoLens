import argparse
import pandas as pd
from Bio import Entrez
import time

def get_species_synonyms(species_name, email, max_retries=3):
    """Fetch synonyms and included names of a species from NCBI with retry mechanism."""
    Entrez.email = email  
    retries = 0
    
    while retries < max_retries:
        try:
            handle = Entrez.esearch(db="taxonomy", term=species_name)
            record = Entrez.read(handle)
            handle.close()
            
            if not record["IdList"]:
                print(f"No tax ID found for {species_name}")
                return species_name, species_name, []
            
            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            tax_record = Entrez.read(handle)[0]
            handle.close()
            
            accepted_name = tax_record.get("ScientificName", species_name)
            synonyms = tax_record.get("OtherNames", {}).get("Synonym", [])
            includes = tax_record.get("OtherNames", {}).get("Includes", [])
            
            all_synonyms = synonyms + includes  
            return species_name, accepted_name, all_synonyms
        
        except Exception as e:
            retries += 1
            print(f"Error fetching data for {species_name} (attempt {retries}/{max_retries}): {e}")
            time.sleep(2 ** retries)  
    
    print(f"Failed to fetch data for {species_name} after {max_retries} attempts.")
    return species_name, species_name, []

def process_species_file(input_file, output_file, email):
    with open(input_file, 'r', encoding='utf-8') as file:
        species_list = [line.strip() for line in file if line.strip()]
    
    results = []
    max_synonyms = 0
    
    for species in species_list:
        time.sleep(0.5)  # Evitar sobrecarga de la API
        eid_sp, accepted_name, synonyms = get_species_synonyms(species, email)
        eid_value = eid_sp if eid_sp != accepted_name else ""
        row = [accepted_name, eid_value] + synonyms
        max_synonyms = max(max_synonyms, len(synonyms))
        results.append(row)
    
    headers = ["Accepted Scientific Name", "EID_sp"] + [f"Synonym {i+1}" for i in range(max_synonyms)]
    
    for row in results:
        row.extend([""] * (max_synonyms - (len(row) - 2)))
    
    df = pd.DataFrame(results, columns=headers)
    df.to_excel(output_file, index=False)
    print(f"File saved: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch species synonyms from NCBI and save to Excel.")
    parser.add_argument("input_file", help="Path to the input text file containing species names.")
    parser.add_argument("output_file", help="Path to the output Excel file.")
    parser.add_argument("email", help="Email address required by NCBI Entrez API.")
    
    args = parser.parse_args()
    process_species_file(args.input_file, args.output_file, args.email)
