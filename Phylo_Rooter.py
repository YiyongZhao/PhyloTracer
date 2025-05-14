from Ortho_Retriever import *
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

import numpy as np
from scipy.stats import norm

def calculate_RF_distance(Phylo_t_OG_L: object, sptree: object) -> int:
    """
    Calculate the Robinson-Foulds (RF) distance between a gene tree and a species tree.

    Args:
        Phylo_t_OG_L (object): The gene tree object.
        sptree (object): The species tree object.

    Returns:
        int: The RF distance.
    """
    for leaf in Phylo_t_OG_L:
        leaf.name = leaf.name.split('_')[0]
    RF = Phylo_t_OG_L.robinson_foulds(sptree)[0]
    return RF

def has_multiple_species(node: object) -> bool:
    """
    Determine whether a node contains genes from multiple species.

    Args:
        node (object): The tree node object, should support calculate_species_num.

    Returns:
        bool: True if the node contains genes from more than one species, otherwise False.
    """
    return calculate_species_num(node) > 1

def get_traversal_string(node: object) -> str:
    """
    Recursively traverse the node and generate a string of node names separated by '-'.

    Args:
        node (object): The tree node object.

    Returns:
        str: A string of node names separated by '-'.
    """
    if has_multiple_species(node):
        return f"{node.name}-{get_traversal_string(node.children[0])}-{get_traversal_string(node.children[1])}"
    return node.name

def generate_rerooted_trees(tree: object, node_name_list: list) -> list:
    """
    Generate a list of rerooted tree objects by setting each node in node_name_list as the outgroup.

    Args:
        tree (object): The original tree object, should support copy, set_outgroup, and name attributes.
        node_name_list (list): List of node names to be used as outgroups.

    Returns:
        list: A list of rerooted tree objects.
    """
    rerooted_trees = []
    for node_name in node_name_list:
        tree_copy = tree.copy('newick')
        outgroup_node = tree_copy & node_name
        tree_copy.set_outgroup(outgroup_node)
        tree_copy.name = node_name
        rerooted_trees.append(tree_copy)
    return rerooted_trees

def get_all_rerooted_trees(tree: object) -> list:
    """
    Get a list of rerooted tree objects after rooting them at each eligible node.

    Args:
        tree (object): The input tree object.

    Returns:
        list: A list of rerooted tree objects.
    """
    tree = num_tre_node(tree)
    node_name_list = get_traversal_string(tree).split('-')
    # Remove root node names from the list
    node_name_list = [name for name in node_name_list if not (tree & name).is_root()]
    rerooted_trees = generate_rerooted_trees(tree, node_name_list)
    return rerooted_trees

def rename_output_tre(tree: object, name_mapping: dict, tree_id: str, output_dir: str) -> None:
    """
    Restore the original gene names in the phylogenetic tree and save it to a file.

    Args:
        tree (object): The phylogenetic tree object.
        name_mapping (dict): Mapping from new gene names to original gene names.
        tree_id (str): Identifier for the tree.
        output_dir (str): Directory to save the output tree file.

    Returns:
        object: The updated phylogenetic tree object.
    """
    for node in tree.traverse():
        if node.name in name_mapping:
            node.name = name_mapping[node.name]

    tree_str=tree.write(format=0)
    write_tree_to_newick(tree_str,tree_id,output_dir)
    
def calculate_tree_statistics(
    tree: object,
    voucher_to_taxa: dict,
    gene_to_new_name: dict,
    tree_id: str,
    tree_path: str,
    renamed_length_dict: dict,
    new_name_to_gene: dict,
    renamed_species_tree: object
) -> tuple:
    """
    Calculate various statistics for a given phylogenetic tree.

    Args:
        tree (object): The phylogenetic tree object.
        voucher_to_taxa (dict): Mapping from voucher to taxa.
        gene_to_new_name (dict): Mapping from gene to new gene name.
        tree_id (str): Identifier for the tree.
        tree_path (str): Path to the tree file.
        renamed_length_dict (dict): Mapping for renamed branch lengths.
        new_name_to_gene (dict): Mapping from new gene name to original gene.
        renamed_species_tree (object): The renamed species tree object.

    Returns:
        tuple: (deep, var, RF, GD, species_overlap)
    """
    up_clade = tree.children[1]
    down_clade = tree.children[0]
    if len(up_clade.get_leaf_names()) > len(down_clade.get_leaf_names()):
        var = abs(compute_tip_to_root_branch_length_variance(up_clade) - compute_tip_to_root_branch_length_variance(down_clade))
        deep = get_max_deepth(down_clade,renamed_species_tree)
    else:
        var = abs(compute_tip_to_root_branch_length_variance(down_clade) - compute_tip_to_root_branch_length_variance(up_clade))
        deep = get_max_deepth(up_clade,renamed_species_tree)
    if len(get_species_list(tree)) == len(get_species_set(tree)):
        RF = calculate_RF_distance(tree, renamed_species_tree)
    else:
        principal_gene_set, filtered_offcut_ev_seqs = offcut_tre(tree, renamed_length_dict)
        minor_orthologs = []
        minor_orthologs = iterator(filtered_offcut_ev_seqs, tree, gene_to_new_name, minor_orthologs, tree_path, renamed_length_dict)
        ordered_name_OG_list = rename_OGs_tre_name(principal_gene_set, minor_orthologs, tree_id)
        RF = 0
        for tre_name, OG_set in ordered_name_OG_list:
            OG_list = [new_name_to_gene[OG] for OG in OG_set]
            phylo_tree_0 = read_phylo_tree(tree_path)
            phylo_tree = root_tre_with_midpoint_outgroup(phylo_tree_0)
            phylo_tree_OG_list = extract_tree(OG_list, phylo_tree)
            phylo_tree_OG_list = rename_input_tre(phylo_tree_OG_list, gene_to_new_name)
            RF += calculate_RF_distance(phylo_tree_OG_list, renamed_species_tree)

    tre_ParaL, GF_leaves_S = find_tre_dup(tree)
    GD = len(tre_ParaL) if tre_ParaL else 0
    species_overlap = calculate_species_overlap(tree)
    return deep, var, RF, GD, species_overlap


def root_main(
    tree_dict: dict,
    gene_to_new_name: dict,
    renamed_length_dict: dict,
    new_name_to_gene: dict,
    renamed_species_tree: object,
    voucher_to_taxa: dict
) -> None:
    """
    Main function for rerooting gene trees, calculating statistics, and selecting the best tree.

    Args:
        tree_dict (dict): Mapping from tree ID to tree file path.
        gene_to_new_name (dict): Mapping from gene to new gene name.
        renamed_length_dict (dict): Mapping for renamed branch lengths.
        new_name_to_gene (dict): Mapping from new gene name to original gene.
        renamed_species_tree (object): The renamed species tree object.
        voucher_to_taxa (dict): Mapping from voucher to taxa.
    """
    dir_path = os.path.join(os.getcwd(), "rooted_trees/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    pbar = tqdm(total=len(tree_dict), desc="Processing trees", unit="tree")


    stat_matrix = []
    weights = {
        "deep": 0.3, #0.5
        "var": 0.05, #0.1
        "RF": 0.5, #0.4
        "GD": 0.1,
        "species_overlap": 0.05
    }

    weights1 = {
        "deep": 0.5, #0.5
        "var": 0.1, #0.1
        "RF": 0.4, #0.4
    }

    try:
        for tree_id, tree_path in tree_dict.items():
            pbar.set_description(f"Processing {tree_id}")
            tree_stats = []
            tree_objects = {}
            
            Phylo_t0 = read_phylo_tree(tree_path)
            Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
            Phylo_t2 = rename_input_tre(Phylo_t1, gene_to_new_name)

            if len(get_species_set(Phylo_t2)) ==1 or len(get_species_list(Phylo_t2))<=3:
                tree_str = Phylo_t2.write()
                rename_output_tre(Phylo_t2, new_name_to_gene, tree_id, dir_path)
                pbar.update(1)
                continue

            root_list = get_all_rerooted_trees(Phylo_t2)
            if root_list:
                root_list = root_list[1:]  # Remove the first (original) root

            for n, tree in enumerate(root_list):
                t1 = rename_input_tre(tree, new_name_to_gene)
                t1.ladderize()
                t1.sort_descendants("support")

                tree_key = f"{tree_id}_{n+1}"
                tree_objects[tree_key] = t1
                
                deep, var, RF, GD, species_overlap = calculate_tree_statistics(
                    tree, voucher_to_taxa, gene_to_new_name, tree_id, tree_path,
                    renamed_length_dict, new_name_to_gene, renamed_species_tree
                )
                
                tree_stats.append({
                    "Tree": tree_key,
                    "deep": deep,
                    "var": var,
                    "RF": RF,
                    "GD": GD,
                    "species_overlap": species_overlap
                })
            
            current_df = pd.DataFrame(tree_stats)

            is_multi_copy = len(get_species_list(Phylo_t2)) != len(get_species_set(Phylo_t2))

            used_weights = weights if is_multi_copy else weights1


            def normalize_and_score(df, weights):
                norm_deep = (df["deep"] - df["deep"].min()) / (df["deep"].max() - df["deep"].min()) if df["deep"].max() != df["deep"].min() else 0
                norm_var = (df["var"] - df["var"].min()) / (df["var"].max() - df["var"].min()) if df["var"].max() != df["var"].min() else 0
                norm_RF = (df["RF"] - df["RF"].min()) / (df["RF"].max() - df["RF"].min()) if df["RF"].max() != df["RF"].min() else 0
                norm_GD = (df["GD"] - df["GD"].min()) / (df["GD"].max() - df["GD"].min()) if df["GD"].max() != df["GD"].min() else 0
                norm_dup_species_overlap = (df["species_overlap"] - df["species_overlap"].min()) / (df["species_overlap"].max() - df["species_overlap"].min()) if df["species_overlap"].max() != df["species_overlap"].min() else 0
                df["weighted_norm_deep"] = norm_deep * weights.get("deep", 0)
                df["weighted_norm_var"] = norm_var * weights.get("var", 0)
                df["weighted_norm_RF"] = norm_RF * weights.get("RF", 0)
                df["weighted_norm_GD"] = norm_GD * weights.get("GD", 0)
                df["weighted_norm_dup_species_overlap"] = norm_dup_species_overlap * weights.get("species_overlap", 0)
                score = (
                    df["weighted_norm_deep"] +
                    df["weighted_norm_var"] +
                    df["weighted_norm_RF"] +
                    df["weighted_norm_GD"] -
                    df["weighted_norm_dup_species_overlap"]
                )
                return score

            current_df["score"] = normalize_and_score(current_df, used_weights)

            for i, stat in enumerate(tree_stats):
                stat["weighted_norm_deep"] = current_df.iloc[i].get("weighted_norm_deep", None)
                stat["weighted_norm_var"] = current_df.iloc[i].get("weighted_norm_var", None)
                stat["weighted_norm_RF"] = current_df.iloc[i].get("weighted_norm_RF", None)
                stat["weighted_norm_GD"] = current_df.iloc[i].get("weighted_norm_GD", None)
                stat["weighted_norm_species_overlap"] = current_df.iloc[i].get("weighted_norm_species_overlap", None)
                stat["score"] = current_df.iloc[i]["score"]

            best_tree_row = current_df.loc[current_df["score"].idxmin()]
            best_tree_key = best_tree_row["Tree"]
            best_tree = tree_objects[best_tree_key]
            best_tree_str = best_tree.write(format=0)
            write_tree_to_newick(best_tree_str, f"{tree_id}", dir_path)

            stat_matrix.extend(tree_stats)
            pbar.update(1)

        stat_df = pd.DataFrame(stat_matrix)
        if not stat_df.empty:
            stat_df["tree_id"] = stat_df["Tree"].apply(lambda x: "_".join(x.split("_")[:-1]))
            stat_df = stat_df.sort_values(by=["tree_id", "score"], ascending=[True, True])
            stat_df.to_csv("stat_matrix3.csv", index=False)
            
    finally:
        pbar.close()
