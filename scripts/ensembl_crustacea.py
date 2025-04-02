#!/usr/bin/env python3

"""
Author: Claudia Restrepo-Ortiz
Date: 2023-02

Usage description:
This script checks the taxonomic classification of arthropod hosts
by querying the Ensembl API.

It requires a CSV input file containing at least two columns:
 - 'Carrier classification': must include information about the host classification.
 - 'Carrier': must contain the species name.
 
For each species classified as an arthropod, the script queries the Ensembl API
to determine if it belongs to the class 'Crustacea', generating a new column
'Crustacean' that will indicate 'yes', 'no', or 'not found' based on the result.
Finally, the results are saved to an Excel file.

python3 script_name.py input_file.csv output_file.xlsx


"""


import requests
import pandas as pd
import argparse
import time

# Function to retrieve taxonomic classification from Ensembl API with retries
def search_tax(sp_name, max_retries=10, delay=2):
    server = "https://rest.ensembl.org"
    ext = "/taxonomy/classification/" + sp_name + "?"
    retries = 0
    
    while retries < max_retries:
        try:
            # Make the request to Ensembl REST API
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})
            
            if r.ok:  # Check if the request was successful
                taxa_dico = r.json()
                SN = "scientific_name"
                TAX = "Crustacea"

                # Check if Crustacea is in the taxonomy list
                for d in taxa_dico:
                    if d.get(SN) == TAX:
                        # Return 'yes' if the species is classified as crustacean
                        return 'yes'
                # Return 'no' if Crustacea is not found in the taxonomy list
                return 'no'
            else:
                print(f"Species '{sp_name}' not found, retrying... ({retries+1}/{max_retries})")
                retries += 1
                time.sleep(delay)  
        except Exception as e:
            print(f"Error occurred: {e}, retrying... ({retries+1}/{max_retries})")
            retries += 1
            time.sleep(delay)
    
    return 'not found'

def main(input_csv, output_file):
    df = pd.read_csv(input_csv)

    df['Crustacean'] = 'no'

    for i in range(len(df)):
        if df['Carrier classification'][i] == 'Arthropod':  # Check if 'Carrier classification' is 'Arthropod'
            sp = df['Carrier'][i]  # Use 'Carrier' for the species name
            crustacean = search_tax(sp)  
            df.at[i, 'Crustacean'] = crustacean  # Update 'Crustacean' column

    df.to_excel(output_file, index=False)
    print(f"Processed data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check taxonomic classification of arthropod hosts.")
    parser.add_argument("input_csv", help="Path to the input CSV file")
    parser.add_argument("output_file", help="Path to the output Excel file")
    
    args = parser.parse_args()
    
    main(args.input_csv, args.output_file)
