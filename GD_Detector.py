from __init__ import *
from collections import Counter

def calculate_depth(node_a: object, node_b: object) -> int:
    """
    Calculate the topological distance between two nodes in a phylogenetic tree.
    Args:
        node_a (object): The first node.
        node_b (object): The second node.
    Returns:
        int: The topological distance between node_a and node_b.
    """
    if node_a in node_b.iter_ancestors() or node_b in node_a.iter_ancestors():
        distance = node_a.get_distance(node_b, topology_only=True) + 2
        return distance
    common_ancestor = node_a.get_common_ancestor(node_b)
    return abs(node_a.get_distance(common_ancestor, topology_only=True) - \
               node_b.get_distance(common_ancestor, topology_only=True))

def map_gene_tree_to_species(species_set: set, species_tree: object) -> object:
    """
    Map a set of species from the gene tree to the corresponding node in the species tree.
    Args:
        species_set (set): Set of species names from the gene tree.
        species_tree (object): The species tree object.
    Returns:
        object: The corresponding node in the species tree.
    """
    if len(species_set) != 1:
        clade = species_tree.get_common_ancestor(species_set)
    else:
        clade = species_tree & list(species_set)[0]
    return clade


def are_sister_supports_greater_than_num(sister_a:object, sister_b:object, min_support: float = 50)->bool:
    """
    Determine whether both sister nodes have support values greater than a given threshold.

    Args:
        sister_a (object): The first sister node, should have a 'support' attribute.
        sister_b (object): The second sister node, should have a 'support' attribute.
        min_support (float or int): The minimum support threshold.

    Returns:
        bool: True if both nodes' support values are greater than min_support, otherwise False.
    """
    support_a = getattr(sister_a, 'support', 0)
    support_b = getattr(sister_b, 'support', 0)
    return support_a > min_support and support_b > min_support

def find_dup_node(
    gene_tree: object,
    species_tree: object,
    gd_support: int = 50,
    clade_support: int = 0,
    dup_species_num: int = 2,
    dup_species_percent: int = 0,
    max_topology_distance: int = 1
) -> list:
    """
    Find duplication nodes in a gene tree based on evolutionary events and various filtering criteria.

    Args:
        gene_tree (object): The gene tree object to analyze.
        species_tree (object): The reference species tree object.
        gd_support (int): Minimum support value for a duplication node to be considered (default: 50).
        clade_support (int): Minimum support value for sister clades (default: 0).
        dup_species_num (int): Minimum number of duplicated species required (default: 2).
        dup_species_percent (int): Minimum percentage of duplicated species required (default: 0).
        max_topology_distance (int): Maximum allowed topological distance between mapped child nodes in the species tree (default: 1).

    Returns:
        list: A list of duplication node objects that meet all criteria.
    """
    dup_node_list = []
    events = gene_tree.get_descendant_evol_events()
    for event in events:
        if event.etype == "D":
            node_names = ",".join(event.in_seqs) + ',' + ",".join(event.out_seqs)
            event_node_name_list = node_names.split(',')
            common_ancestor_node = gene_tree.get_common_ancestor(event_node_name_list)
            child_a, child_b = common_ancestor_node.get_children()
            species_set = get_species_set(common_ancestor_node)
            mapped_species_node = map_gene_tree_to_species(species_set, species_tree)
            common_ancestor_node.add_feature('map', mapped_species_node.name)
            
            # 检查重复节点支持值
            if judge_support(common_ancestor_node.support, gd_support):
                child_a, child_b = common_ancestor_node.get_children()
                
                # 检查子分支支持值
                if child_a.support >= clade_support and child_b.support >= clade_support:
                    mapped_a = map_gene_tree_to_species(get_species_set(child_a), species_tree)
                    mapped_b = map_gene_tree_to_species(get_species_set(child_b), species_tree)
                    
                    if len(get_species_set(common_ancestor_node))==1:
                        dup_node_list.append(common_ancestor_node)
                    else:
                        # 计算重复物种数量和百分比
                        dup_sps = count_common_elements(get_species_set(child_a), get_species_set(child_b))
                        dup_percent = dup_sps / len(get_species_set(common_ancestor_node))
                        # 检查重复物种数量和百分比是否满足条件
                        if dup_sps >= dup_species_num and dup_percent >= dup_species_percent:
                            # 检查拓扑距离
                            if species_tree.get_distance(mapped_a, mapped_b, topology_only=True) <= max_topology_distance:
                                dup_node_list.append(common_ancestor_node)
    return dup_node_list

def write_gd_result(
    output_file: str,
    tree_dict: dict,
    gd_support: int,
    clade_support: int,
    dup_species_percent: float,
    dup_species_num: int,
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
        tree_dict (dict): Dictionary mapping tree IDs to tree file paths.
        gd_support (int): Minimum support value for a duplication node to be considered.
        clade_support (int): Minimum support value for sister clades.
        dup_species_percent (float): Minimum percentage of duplicated species required.
        dup_species_num (int): Minimum number of duplicated species required.
        species_tree (object): The reference species tree object.
        gene_to_new_name (dict): Mapping from gene names to new names.
        new_name_to_gene (dict): Mapping from new names back to gene names.
        voucher_to_taxa (dict): Mapping from voucher IDs to taxa.
        max_topology_distance (int): Maximum allowed topological distance between mapped child nodes in the species tree.

    Returns:
        None
    """
    print(f"=== Gene Duplication Detection Configuration ===")
    print(f"GD node support threshold: {gd_support}")
    print(f"Subclade support threshold: {clade_support}")
    print(f"Minimum duplicated species count: {dup_species_num}")
    print(f"Duplicated species percentage threshold: {dup_species_percent}")
    print(f"Maximum topology distance: {max_topology_distance}")
    print(f"===============================================")
    with open(output_file, 'w') as file:
        file.write('#tree_ID\tgd_id\tgd_support\tgene1\tgene2\tlevel\tspecies\tGD_dup_sps\tdup_ratio\tgd_type\tcomment\n')
        gd_num = 1
        gd_type_dict = {}
        for tree_id, tree_path in tree_dict.items():
            gene_tree = read_phylo_tree(tree_path)
            gene_tree = rename_input_tre(gene_tree, gene_to_new_name)
            if len(gene_tree.children)!=2:
                print(f'{tree_id} is not a binary tree')
                continue
            num_tre_node(gene_tree)
            dup_node_list = find_dup_node(gene_tree, species_tree, gd_support, clade_support, dup_species_num, dup_species_percent, max_topology_distance)
            for node in dup_node_list:
                species_set = get_species_set(node)
                mapped_species_node = map_gene_tree_to_species(species_set, species_tree)
                clade = node
                parent = clade.map
                child_a, child_b = clade.get_children()
                dup_sps =  count_common_elements(get_species_set(child_a), get_species_set(child_b))
                dup_percent = dup_sps / len(get_species_set(clade))
                model = get_model(clade, species_tree)
                gene_pairs = gene_pair(clade)
                null = 'null'
                if not mapped_species_node.is_leaf():
                    if mapped_species_node.name in gd_type_dict:
                        gd_type_dict[mapped_species_node.name].append(model)
                    else:
                        gd_type_dict[mapped_species_node.name]=[model]

                
                for gene_pair_str in gene_pairs:
                    file.write(str(tree_id) + '\t' + str(gd_num) + '\t')
                    gene_a, gene_b = gene_pair_str.split('-')
                    if gene_a.split('_')[0] == gene_b.split('_')[0]:
                        species_code = gene_a.split('_')[0]
                    else:
                        if gene_a == 'null':
                            species_code = gene_b.split('_')[0]
                        else:
                            species_code = gene_a.split('_')[0]
                    file.write(str(clade.support) + '\t' + new_name_to_gene.get(gene_a, null) + '\t' + new_name_to_gene.get(gene_b, null) + '\t' + voucher_to_taxa.get(parent, parent) + '\t' + voucher_to_taxa[species_code] + '\t' + str(dup_sps) + '\t' + str(round(dup_percent, 2)) + '\t' + model + '\t-\t\n')
                gd_num += 1
            
        merged_data=count_elements_in_lists(gd_type_dict)
        df = pd.DataFrame.from_dict(merged_data, orient='index').fillna(0)

        df.to_csv('gd_type.csv')

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
    result = merge_and_filter_types(counted_data, merge_map)
    return result

def get_model(clade: object, sptree: object) -> str:
    """
    Assign labels to species clades and generate a model string based on the gene clade and species tree.

    Args:
        clade (object): The gene clade object, whose leaves contain species information in their names.
        sptree (object): The species tree object, supporting methods like get_common_ancestor and get_children.

    Returns:
        str: A string representing the model type, e.g., 'AAB', 'ABB', etc.
    """
    sps = get_species_list(clade)
    sps_clade = sptree.get_common_ancestor(set(sps))

    sps_clade_a = sps_clade.get_children()[0] if sps_clade.get_children() else None
    sps_clade_b = sps_clade.get_children()[1] if sps_clade.get_children() and len(sps_clade.get_children()) > 1 else None

    if sps_clade_a.is_leaf():
        sps_clade_a.add_feature('label', 'Aa')
    else:
        sps_clade_a_1 = sps_clade_a.get_children()[0] if sps_clade_a.get_children() else None
        sps_clade_a_2 = sps_clade_a.get_children()[1] if sps_clade_a.get_children() and len(sps_clade_a.get_children()) > 1 else None
        for leaf in sps_clade_a_1:
            leaf.add_feature('label', 'A')
        for leaf in sps_clade_a_2:
            leaf.add_feature('label', 'a')

    if sps_clade_b.is_leaf():
        sps_clade_b.add_feature('label', 'Bb')
    else:
        sps_clade_b_1 = sps_clade_b.get_children()[0] if sps_clade_b.get_children() else None
        sps_clade_b_2 = sps_clade_b.get_children()[1] if sps_clade_b.get_children() and len(sps_clade_b.get_children()) > 1 else None
        for leaf in sps_clade_b_1:
            leaf.add_feature('label', 'B')
        for leaf in sps_clade_b_2:
            leaf.add_feature('label', 'b')

    for node in clade:
        species = node.name.split('_')[0]
        clade1 = sps_clade & species
        if clade1:
            node.add_feature('label', clade1.label)

    up_clade = ''
    for node in clade.get_children()[0]:
        up_clade += node.label
    up_clade = up_clade + '<=>'
    for node in clade.get_children()[1]:
        up_clade += node.label
    clade_up = set(up_clade.split('<=>')[0])
    clade_down = set(up_clade.split('<=>')[1])
    clade_up_1 = ''.join(clade_up)
    clade_up_1_1 = ''.join(sorted(clade_up_1, key=lambda x: (x.lower(), x.isupper())))

    clade_down_1 = ''.join(clade_down)
    clade_down_1_1 = ''.join(sorted(clade_down_1, key=lambda x: (x.lower(), x.isupper())))
    clade_model = clade_up_1_1 + '<=>' + clade_down_1_1

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
    while i < len(input_string):
        if i < len(input_string) - 1 and ((input_string[i] == 'A' and input_string[i+1] == 'a') or (input_string[i] == 'a' and input_string[i+1] == 'A')):
            result.append('A')
            i += 2
        elif i < len(input_string) - 1 and ((input_string[i] == 'B' and input_string[i+1] == 'b') or (input_string[i] == 'b' and input_string[i+1] == 'B')):
            result.append('B')
            i += 2
        else:
            if input_string[i] in ['A', 'B', 'a', 'b']:
                result.append('X')
            else:
                result.append(input_string[i])
            i += 1

    return ''.join(result)
    
def get_model_dic(interspecies_node_list: list, genetree: object, sptree: object) -> dict:
    """
    Generate a dictionary mapping model types to interspecies nodes based on the gene tree and species tree.

    Args:
        interspecies_node_list (list): A list of interspecies node identifiers.
        genetree (object): The gene tree object, supporting '&' operation to extract clades.
        sptree (object): The species tree object, used to determine the model type for each clade.

    Returns:
        dict: A dictionary where keys are model types (as determined by get_model) and values are lists of interspecies node identifiers.
    """
    model_dic = {}
    for node_id in interspecies_node_list:
        clade = genetree & node_id
        model_type = get_model(clade, sptree)
        model_dic.setdefault(model_type, []).append(node_id)
    return model_dic

def get_empty_count_dict(sptree: object) -> dict:
    """
    Generate an empty count dictionary for each node in the species tree.

    Args:
        sptree (object): The species tree object, each node should have a name attribute.

    Returns:
        dict: A dictionary with node names as keys and zero as values.
    """
    empty_count_dic = {node.name: 0 for node in sptree.traverse()}
    return empty_count_dic

def judge_support(support: float, support_value: float) -> bool:
    """
    Judge whether the support value meets the threshold, supporting both proportion (0-1) and percentage (0-100) formats.

    Args:
        support (float): The support value to be judged, can be a proportion or percentage.
        support_value (float): The threshold value, can be a proportion (0-1) or percentage (0-100).

    Returns:
        bool: True if the support meets or exceeds the threshold, False otherwise.
    """
    if support <= 1 and 0.5 <= support_value <= 1:
        if support >= support_value:
            return True
        else:
            return False
    elif support <= 1 and 50 <= support_value <= 100:
        support_value = support_value / 100
        if support >= support_value:
            return True
        else:
            return False
    elif support > 1 and 0.5 <= support_value <= 1:
        support_value = support_value * 100
        if support >= support_value:
            return True
        else:
            return False
    elif support > 1 and 50 <= support_value <= 100:
        if support >= support_value:
            return True
        else:
            return False

def gene_pair(clade:object) -> set:
    """
    Generate gene pairs from two child clades based on matching species codes in leaf names.

    Args:
        clade: A clade object with a method get_children(), whose children have get_leaf_names().

    Returns:
        set: A set of gene pair strings in the format 'geneA-geneB', 'geneA-null', or 'null-geneB'.
    """
    result_pairs = set()

    children = clade.get_children()
    child1, child2 = sorted(children, key=lambda c: len(c.get_leaf_names()), reverse=True)
    leaves1, leaves2 = child1.get_leaf_names(), child2.get_leaf_names()

    for tip1 in leaves1:
        matching_tips = [tip2 for tip2 in leaves2 if tip1.split('_')[0] == tip2.split('_')[0]]
        if matching_tips:
            result_pairs.update(f"{tip1}-{tip2}" for tip2 in matching_tips)
        else:
            result_pairs.add(f"{tip1}-null")

    for tip2 in leaves2:
        if all(tip2.split('_')[0] != tip1.split('_')[0] for tip1 in leaves1):
            result_pairs.add(f"null-{tip2}")

    return result_pairs

if __name__ == "__main__":
    support=50
    dup_species_percent = 0.5
    dup_species_num = 2
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap.txt")
    sptree=PhyloTree('30sptree.nwk')
    sptree=rename_species_tree(sptree, voucher2taxa_dic)
    num_tre_node(sptree)
    tre_dic=read_and_return_dict('GF.txt')
    filename = 'result.txt'
    sp_dic=[]
    write_gd_result(sp_dic,filename, tre_dic, support,dup_species_percent, dup_species_num,sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic)
   
