import os
import pandas as pd
from patholens.utils import load_species_from_txt, load_excel_data, find_species_in_excel
from patholens.config import INPUT_DIR


def extract_final_species(output_file,group,dir_group):
   
    
    txt_file=os.path.join(dir_group, f"FINAL_unique_taxonomies_{group}.txt")
    excel_file=os.path.join(INPUT_DIR, group, f"{group}_Pathogen_TaxSyn_List.xlsx")

    species_list = load_species_from_txt(txt_file)
    df = load_excel_data(excel_file)
    results = find_species_in_excel(species_list, df)

    output_df = pd.DataFrame(results, columns=['Species Found in DB','Accepted Scientific Name' ])
    grouped_df = output_df.groupby('Accepted Scientific Name')['Species Found in DB'].apply(lambda x: ', '.join(sorted(set(x)))).reset_index()
    grouped_df.to_excel(output_file, index=False)

    print("Total species in the curated Database: ", len(grouped_df))
