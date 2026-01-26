from __init__ import *
from collections import Counter
import pandas as pd

def write_gene_duplication_results(
    output_file: str,
    tree_paths: dict,
    duplication_support_threshold: int,
    subclade_support_threshold: int,
    duplicated_species_percentage_threshold: float,
    duplicated_species_count_threshold: int,
    species_tree: object,
    gene_to_new_name: dict,
    new_name_to_gene: dict,
    voucher_to_taxa: dict,
    max_topology_distance: int
):
    """
    Write gene duplication detection results to a file for a set of gene trees.

    Args:
        output_file (str): Path to the output file.
        tree_paths (dict): Dictionary mapping tree IDs to tree file paths.
        duplication_support_threshold (int): Minimum support value for a duplication node to be considered.
        subclade_support_threshold (int): Minimum support value for sister clades.
        duplicated_species_percentage_threshold (float): Minimum percentage of duplicated species required.
        duplicated_species_count_threshold (int): Minimum number of duplicated species required.
        species_tree (object): The reference species tree object.
        gene_to_new_name (dict): Mapping from gene names to new names.
        new_name_to_gene (dict): Mapping from new names back to gene names.
        voucher_to_taxa (dict): Mapping from voucher IDs to taxa.
        max_topology_distance (int): Maximum allowed topological distance between mapped child nodes in the species tree.

    Returns:
        None
    """
    print("=== Gene Duplication Detection Configuration ===")
    print(f"GD node support threshold: {duplication_support_threshold}")
    print(f"Subclade support threshold: {subclade_support_threshold}")
    print(f"Minimum duplicated species count: {duplicated_species_count_threshold}")
    print(f"Duplicated species percentage threshold: {duplicated_species_percentage_threshold}")
    print(f"Maximum topology distance: {max_topology_distance}")
    print("===============================================")

    gd_type_dict = {}
    gd_num = 1
    null = 'null'

    with open(output_file, 'w') as file:
        file.write('#tree_ID\tgd_id\tgd_support\tgene1\tgene2\tlevel\tspecies\tGD_dup_sps\tdup_ratio\tgd_type\tcomment\n')

        for tree_id, tree_path in tree_paths.items():
            gene_tree = read_phylo_tree(tree_path)
            gene_tree = rename_input_tre(gene_tree, gene_to_new_name)
            if len(gene_tree.children) != 2:
                print(f'{tree_id} is not a binary tree')
                continue
            num_tre_node(gene_tree)
            dup_node_list = find_dup_node(
                gene_tree,
                species_tree,
                duplication_support_threshold,
                subclade_support_threshold,
                duplicated_species_count_threshold,
                duplicated_species_percentage_threshold,
                max_topology_distance
            )
            for node in dup_node_list:
                clade = node
                species_set = get_species_set(clade)
                mapped_species_node = map_gene_tree_to_species(species_set, species_tree)
                child_a, child_b = clade.get_children()
                child_a_species_set = get_species_set(child_a)
                child_b_species_set = get_species_set(child_b)
                dup_species_count = count_common_elements(child_a_species_set, child_b_species_set)
                dup_percentage = dup_species_count / len(species_set) if species_set else 0
                model = get_model(clade, species_tree)
                gene_pairs = gene_pair(clade)
                parent = clade.map if hasattr(clade, 'map') else None

                if not mapped_species_node.is_leaf():
                    if mapped_species_node.name in gd_type_dict:
                        gd_type_dict[mapped_species_node.name].append(model)
                    else:
                        gd_type_dict[mapped_species_node.name] = [model]

                for gene_pair_str in gene_pairs:
                    gene_a, gene_b = gene_pair_str.split('-')
                    if gene_a.split('_')[0] == gene_b.split('_')[0]:
                        species_code = gene_a.split('_')[0]
                    else:
                        species_code = gene_b.split('_')[0] if gene_a == 'null' else gene_a.split('_')[0]
                    level = voucher_to_taxa.get(parent, parent) if parent else 'unknown'
                    species = voucher_to_taxa.get(species_code, species_code)
                    file.write(
                        f"{tree_id}\t{gd_num}\t{clade.support}\t"
                        f"{new_name_to_gene.get(gene_a, null)}\t{new_name_to_gene.get(gene_b, null)}\t"
                        f"{level}\t{species}\t{dup_species_count}\t{round(dup_percentage, 2)}\t"
                        f"{model}\t-\n"
                    )
                gd_num += 1

    merged_data = count_elements_in_lists(gd_type_dict)
    df = pd.DataFrame.from_dict(merged_data, orient='index').fillna(0)
    df.to_csv('gd_type.tsv', sep='\t')

def merge_and_filter_types(
    data: dict,
    merge_map: dict
) -> dict:
    """
    Merge and filter event types in the input data according to the merge_map.

    Args:
        data (dict): A dictionary where keys are identifiers and values are Counter objects of event types.
        merge_map (dict): A mapping from new event types to lists of event types to be merged.

    Returns:
        dict: A dictionary with the same keys as data, but with merged event type counts.
    """
    merged_data = {}
    for key, counter in data.items():
        new_counter = Counter()
        for event_type, count in counter.items():
            merged = False
            for new_type, types_to_merge in merge_map.items():
                if event_type in types_to_merge:
                    new_counter[new_type] += count
                    merged = True
                    break
            if not merged and event_type == 'AB<=>AB':
                new_counter['ABAB'] += count
        merged_data[key] = new_counter
    return merged_data

def count_elements_in_lists(data: dict) -> dict:
    """
    Count and merge event types in nested lists using merge_and_filter_types.

    Args:
        data (dict): A dictionary where keys are identifiers and values are lists of event types.

    Returns:
        dict: A dictionary with merged event type counts for each identifier.
    """
    counted_data = {key: Counter(value) for key, value in data.items()}
    merge_map = {
        'ABB': ['AB<=>B', 'B<=>AB', 'XB<=>AB', 'AB<=>XB'],
        'AAB': ['AB<=>A', 'A<=>AB', 'AX<=>AB', 'AB<=>AX']
    }
    return merge_and_filter_types(counted_data, merge_map)

def get_model(clade: object, species_tree: object) -> str:
    """
    Assign labels to species clades and generate a model string based on the gene clade and species tree.

    Args:
        clade (object): The gene clade object, whose leaves contain species information in their names.
        species_tree (object): The species tree object, supporting methods like get_common_ancestor and get_children.

    Returns:
        str: A string representing the model type, e.g., 'AAB', 'ABB', etc.
    """
    species_list = get_species_list(clade)
    species_clade = species_tree.get_common_ancestor(set(species_list))

    if not species_clade.get_children():
        return ''

    species_clade_a, species_clade_b = species_clade.get_children()[:2]

    if species_clade_a.is_leaf():
        species_clade_a.add_feature('label', 'Aa')
    else:
        children = species_clade_a.get_children()[:2]
        if children:
            species_clade_a_1, species_clade_a_2 = children
            for leaf in species_clade_a_1.get_leaves():
                leaf.add_feature('label', 'A')
            for leaf in species_clade_a_2.get_leaves():
                leaf.add_feature('label', 'a')

    if species_clade_b.is_leaf():
        species_clade_b.add_feature('label', 'Bb')
    else:
        children = species_clade_b.get_children()[:2]
        if children:
            species_clade_b_1, species_clade_b_2 = children
            for leaf in species_clade_b_1.get_leaves():
                leaf.add_feature('label', 'B')
            for leaf in species_clade_b_2.get_leaves():
                leaf.add_feature('label', 'b')

    for node in clade.get_leaves():
        species = node.name.split('_')[0]
        species_leaf = species_clade & species
        if species_leaf:
            node.add_feature('label', species_leaf.label)

    children = clade.get_children()
    if len(children) != 2:
        return ''

    up_clade = ''
    up_clade += ''.join(node.label for node in children[0].get_leaves() if hasattr(node, 'label'))
    up_clade += '<=>'
    up_clade += ''.join(node.label for node in children[1].get_leaves() if hasattr(node, 'label'))

    clade_up = set(up_clade.split('<=>')[0])
    clade_down = set(up_clade.split('<=>')[1])
    clade_up_sorted = ''.join(sorted(''.join(clade_up), key=lambda x: (x.lower(), x.isupper())))
    clade_down_sorted = ''.join(sorted(''.join(clade_down), key=lambda x: (x.lower(), x.isupper())))
    clade_model = clade_up_sorted + '<=>' + clade_down_sorted

    return process_string(clade_model)

def process_string(input_string: str) -> str:
    """
    Process the input string by replacing specific character combinations.

    Args:
        input_string (str): The input string to be processed.

    Returns:
        str: The processed string with specific character combinations replaced.
    """
    result = []
    i = 0
    length = len(input_string)
    while i < length:
        if i < length - 1 and input_string[i].upper() == input_string[i+1].upper() and input_string[i].lower() == input_string[i+1].lower():
            result.append(input_string[i].upper())
            i += 2
        else:
            char = input_string[i]
            if char.upper() in ['A', 'B']:
                result.append('X')
            else:
                result.append(char)
            i += 1
    return ''.join(result)

def get_model_dic(interspecies_node_list: list, gene_tree: object, species_tree: object) -> dict:
    """
    Generate a dictionary mapping model types to interspecies nodes based on the gene tree and species tree.

    Args:
        interspecies_node_list (list): A list of interspecies node identifiers.
        gene_tree (object): The gene tree object, supporting '&' operation to extract clades.
        species_tree (object): The species tree object, used to determine the model type for each clade.

    Returns:
        dict: A dictionary where keys are model types (as determined by get_model) and values are lists of interspecies node identifiers.
    """
    model_dict = {}
    for node_id in interspecies_node_list:
        clade = gene_tree & node_id
        if clade:
            model_type = get_model(clade, species_tree)
            model_dict.setdefault(model_type, []).append(node_id)
    return model_dict

def get_empty_count_dict(species_tree: object) -> dict:
    """
    Generate an empty count dictionary for each node in the species tree.

    Args:
        species_tree (object): The species tree object, each node should have a name attribute.

    Returns:
        dict: A dictionary with node names as keys and zero as values.
    """
    return {node.name: 0 for node in species_tree.traverse() if hasattr(node, 'name')}

def gene_pair(clade: object) -> set:
    """
    Generate gene pairs from two child clades based on matching species codes in leaf names.

    Args:
        clade: A clade object with a method get_children(), whose children have get_leaf_names().

    Returns:
        set: A set of gene pair strings in the format 'geneA-geneB', 'geneA-null', or 'null-geneB'.
    """
    result_pairs = set()
    children = clade.get_children()
    if len(children) != 2:
        return result_pairs

    child1, child2 = sorted(children, key=lambda c: len(c.get_leaf_names()), reverse=True)
    leaves1 = child1.get_leaf_names()
    leaves2 = child2.get_leaf_names()

    species_to_genes1 = {}
    for tip in leaves1:
        species = tip.split('_')[0]
        species_to_genes1.setdefault(species, []).append(tip)

    species_to_genes2 = {}
    for tip in leaves2:
        species = tip.split('_')[0]
        species_to_genes2.setdefault(species, []).append(tip)

    all_species = set(species_to_genes1) | set(species_to_genes2)

    for species in all_species:
        genes1 = species_to_genes1.get(species, [])
        genes2 = species_to_genes2.get(species, [])
        if genes1 and genes2:
            for g1 in genes1:
                for g2 in genes2:
                    result_pairs.add(f"{g1}-{g2}")
        elif genes1:
            for g1 in genes1:
                result_pairs.add(f"{g1}-null")
        elif genes2:
            for g2 in genes2:
                result_pairs.add(f"null-{g2}")

    return result_pairs

if __name__ == "__main__":
    duplication_support_threshold = 50
    subclade_support_threshold = 50  # Assuming same as duplication_support_threshold based on original usage
    duplicated_species_percentage_threshold = 0.5
    duplicated_species_count_threshold = 2
    max_topology_distance = 0  # Assuming default value; adjust if needed
    gene_to_new_name, new_name_to_gene, voucher_to_taxa = gene_id_transfer("imap.txt")
    species_tree = PhyloTree('30sptree.nwk')
    species_tree = rename_species_tree(species_tree, voucher_to_taxa)
    num_tre_node(species_tree)
    tree_paths = read_and_return_dict('GF.txt')
    output_file = 'result.txt'
    write_gene_duplication_results(
        output_file,
        tree_paths,
        duplication_support_threshold,
        subclade_support_threshold,
        duplicated_species_percentage_threshold,
        duplicated_species_count_threshold,
        species_tree,
        gene_to_new_name,
        new_name_to_gene,
        voucher_to_taxa,
        max_topology_distance
    )
