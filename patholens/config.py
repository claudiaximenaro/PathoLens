import os

#relative paths to the PathoLens directory
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
INPUT_DIR = os.path.join(DATA_DIR, "input")
OUTPUT_DIR = os.path.join(DATA_DIR, "output")
GENERAL_RESULTS_FILE = os.path.join(OUTPUT_DIR, "all_groups_general_results.csv")


#subdirectory for SILVA database
SILVA_DIR = os.path.join(INPUT_DIR, "SILVA")  # Directory for SILVA FASTA files
os.makedirs(OUTPUT_DIR, exist_ok=True)

# get the output path dynamically for each taxonomic group
def get_group_output_dir(group):
    """
    Generates the output path for a given group.

    Parameters:
        group (str): The name of the group.

    Returns:
        str: Path to the group's output directory.
    """
    group_output_dir = os.path.join(OUTPUT_DIR, group)
    os.makedirs(group_output_dir, exist_ok=True)
    return group_output_dir