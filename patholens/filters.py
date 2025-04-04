from patholens.utils import count_and_extract_taxonomies, write_to_csv_general_results


def compare_taxonomy_genus_species(fasta_file,group):
    discrepancies = []
    discrepancy_count = 0
    unique_taxonomies, complete_taxonomies =count_and_extract_taxonomies(fasta_file)

    for taxonomy in unique_taxonomies:
        parts = taxonomy.split(';')

        if len(parts) < 6:
            discrepancy_count += complete_taxonomies.count(taxonomy)
            discrepancies.append(taxonomy)
        else:
            genus = parts[-2]
            species = parts[-1]
            genus_from_species = species.split()[0]
            if genus != genus_from_species:
                discrepancy_count += complete_taxonomies.count(taxonomy)
                discrepancies.append(taxonomy)
    discrepancies_uniques=set(discrepancies)

    write_to_csv_general_results(discrepancy_count, len(discrepancies_uniques), group, "Genus-Species Discrepancy")


    return discrepancies_uniques, complete_taxonomies


def filter_gen_included(rest_list, group, complete_taxonomies):
    gen_included_output = []
    to_review_output = []
    total_review_count = 0

    def clean_taxonomy(taxonomy):
        return [term for part in taxonomy.split('-') for term in part.split()]

    for line in rest_list:
        columns = line.split(';')
        gen_taxonomy = columns[-2]
        species_full = columns[-1].strip()
        genus = species_full.split()[0]

        normalized_taxonomy = clean_taxonomy(gen_taxonomy)

        if genus in normalized_taxonomy:
            gen_included_output.append(line)
        else:
            to_review_output.append(line)
            total_review_count += complete_taxonomies.count(line)

    write_to_csv_general_results(total_review_count, len(to_review_output), group, "Genus Included")

    return gen_included_output, to_review_output