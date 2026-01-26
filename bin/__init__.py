import logging
import pandas as pd
from ete3 import PhyloTree,Tree,NodeStyle,TreeStyle,TextFace,RectFace
import random
import numpy as np
import os
import shutil
from tqdm import tqdm

def generate_sps_voucher(sps_num: int) -> list:
    """
    Generate a list of unique 3-character vouchers for species labeling.

    Args:
        sps_num (int): Number of unique vouchers to generate.

    Returns:
        list: Sorted list of unique voucher strings.
    """
    characters = [chr(i) for i in range(65, 91)] + [chr(i) for i in range(97, 123)] + [str(i) for i in range(10)]
    unique_strings = set()
    while len(unique_strings) < sps_num:
        unique_strings.add(''.join(random.sample(characters, 3)))
    return sorted(list(unique_strings))

def gene_id_transfer(gene2taxa_list: str) -> dict:
    """
    Transfer gene IDs to new voucher-based IDs and build mapping dictionaries.

    Args:
        gene2taxa_list (str): Path to the gene-to-taxa mapping file.

    Returns:
        tuple: Dictionaries for gene-to-new-name, new-name-to-gene, voucher-to-taxa, and taxa-to-voucher.
    """
    gene2taxa_dic = read_and_return_dict(gene2taxa_list)
    sorted_gene2taxa_dic = dict(sorted(gene2taxa_dic.items(), key=lambda x: x[1]))
    taxa_list = list(set(sorted_gene2taxa_dic.values()))
    taxa2voucher_dic = dict(zip(taxa_list, generate_sps_voucher(len(taxa_list))))
    voucher2taxa_dic = {value: key for key, value in taxa2voucher_dic.items()}
    gene_count = {}
    for species in sorted_gene2taxa_dic.values():
        gene_count[species] = gene_count.get(species, 0) + 1
    new_gene_names = [f"{taxa2voucher_dic[species]}_{i}" for species, count in gene_count.items() for i in range(1, count + 1)]
    gene2new_named_gene_dic = dict(zip(sorted_gene2taxa_dic.keys(), new_gene_names))
    new_named_gene2gene_dic = {value: key for key, value in gene2new_named_gene_dic.items()}
    return gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic
#gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("gene2taxa.list")

def read_and_return_dict(filename: str, separator: str = "\t") -> dict:
    """
    Read a two-column mapping file and return it as a dictionary.

    Args:
        filename (str): Path to the mapping file.
        separator (str): Column separator.

    Returns:
        dict: Mapping from first column to second column.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    try:
        df = pd.read_csv(filename, sep=separator, header=None)
        if df.empty:
            raise ValueError("Mapping file is empty")
        return df.set_index([0])[1].to_dict()
    except Exception as e:
        logging.error(f"Failed to parse mapping file {filename}: {e}")
        raise

def mapp_gene_tree_to_species(sp_set: set, sptree: object) -> object:
    """
    Map a set of species names to their most recent common ancestor (MRCA) in a species tree.

    Args:
        sp_set (set): A set of species names (strings) corresponding to tip names in the species tree.
        sptree (object): The input species tree (ETE3 Tree object).

    Returns:
        object: The tree node representing the most recent common ancestor (MRCA) of all species
                in sp_set. If sp_set contains only one species, the corresponding leaf node is returned.

    Example:
        sp_set = {'Homo_sapiens', 'Mus_musculus'}
        clade_node = mapp_gene_tree_to_species(sp_set, species_tree)
    """
    if len(sp_set) != 1:
        clade = sptree.get_common_ancestor(sp_set)
    else:
        clade = sptree & list(sp_set)[0]
    return clade


def rename_input_tre(Phylo_t: object, gene2new_named_gene_dic: dict) -> object:
    """
    Rename the leaf nodes of a phylogenetic tree according to a mapping dictionary.

    Args:
        Phylo_t (object): The input phylogenetic tree.
        gene2new_named_gene_dic (dict): Mapping from old gene names to new names.

    Returns:
        object: A copy of the tree with renamed nodes.
    """
    Phylo_t1=Phylo_t.copy()
    for node in Phylo_t1.traverse():
        if node.name in gene2new_named_gene_dic:
            node.name = gene2new_named_gene_dic[node.name]
    return Phylo_t1

def read_tree(tree_path: str) -> object:
    """
    Read a tree from a file and return a Tree object.

    Args:
        tree_path (str): Path to the tree file.

    Returns:
        object: Tree object.
    """
    try:
        return Tree(tree_path, format=0)
    except:
        return Tree(tree_path, format=1)

def read_phylo_tree(tree_path: str) -> object:
    """
    Read a phylogenetic tree from a file and return a PhyloTree object.

    Args:
        tree_path (str): Path to the phylogenetic tree file.

    Returns:
        object: PhyloTree object.
    """
    try:
        return PhyloTree(tree_path, format=0)
    except:
        try:
            return PhyloTree(tree_path, format=1)
        except Exception as e:
            print(f"Error reading tree: {tree_path}")
            raise e

def root_tre_with_midpoint_outgroup(Phylo_t: object) -> object:
    """Rooting the phylogenetic tree using the midpoint outgroup method."""
    Phylo_t1 = Phylo_t.copy('newick')
    
    # 检查树是否已经有根
    if is_rooted(Phylo_t1):
        return Phylo_t1
    
    # 检查树的节点数量，少于3个叶节点无需定根
    leaves = Phylo_t1.get_leaves()
    if len(leaves) < 3:
        return Phylo_t1
    
    try:
        mid_node = Phylo_t1.get_midpoint_outgroup()
        
        # 检查中点外群是否是根节点
        if mid_node.is_root():
            # 如果中点外群是根节点，使用第一个叶节点作为外群
            Phylo_t1.set_outgroup(leaves[0])
        else:
            Phylo_t1.set_outgroup(mid_node)
            
    except Exception as e:
        # 如果中点定根失败，使用第一个叶节点作为外群
        print(f"Warning: Midpoint rooting failed ({e}), using first leaf as outgroup")
        if leaves:
            Phylo_t1.set_outgroup(leaves[0])
    
    return Phylo_t1


def is_rooted(Phylo_t: object) -> bool:
    """
    Determine whether the input phylogenetic tree is rooted.

    Args:
        Phylo_t (object): The phylogenetic tree object.

    Returns:
        bool: True if the tree is rooted (i.e., has exactly two children at the root), otherwise False.
    """
    return len(Phylo_t.get_children()) == 2


def num_tre_node(Phylo_t: object) -> object:
    """
    Number the internal nodes of a phylogenetic tree in postorder traversal.

    Args:
        Phylo_t (object): The phylogenetic tree object.

    Returns:
        object: The tree with internal nodes numbered as "N0", "N1", etc.
    """
    i = 0
    for node in Phylo_t.traverse('postorder'):
        if not node.is_leaf():
            node.name = "N" + str(i)
            i += 1
    return Phylo_t

def get_species_list(node: object) -> list:
    """
    Get the list of species names under a given tree node.

    Args:
        node (object): Tree node.

    Returns:
        list: List of species names.
    """
    if node is None:
        return []
    return [leaf.name.split('_')[0] for leaf in node.iter_leaves()]


def find_dup_node_simple(Phylo_t: Tree) -> list:
    """
    Find all duplication nodes in a phylogenetic tree.

    Args:
        Phylo_t (Tree): Phylogenetic tree object.

    Returns:
        list: List of duplication nodes.
    """
    events = Phylo_t.get_descendant_evol_events()
    dup_node_list = []
    for ev in events:
        if ev.etype == "D":
            try:
                in_nodes = [Phylo_t & seq for seq in ev.in_seqs]
                out_nodes = [Phylo_t & seq for seq in ev.out_seqs]
                lca = Phylo_t.get_common_ancestor(in_nodes + out_nodes)
                dup_node_list.append(lca)
            except Exception as e:
                logging.warning(f"Skipping event due to error: {e}")
                continue
    return dup_node_list

def get_species_set(Phylo_t: object) -> set:
    """
    Get the set of unique species names in a phylogenetic tree.

    Args:
        Phylo_t (object): Phylogenetic tree object.

    Returns:
        set: Set of unique species names.
    """
    return set(get_species_list(Phylo_t))

def get_max_deepth(root:object)->int:
    """
    Calculate the maximum depth of a tree.
    Args:
        root (object): The root node of the tree, which should have a 'children' attribute.
    Returns:
        int: The maximum depth of the tree (root counts as level 1).
    """
    if not root:
        return 0
    
    max_child_depth = 0
    for child in root.children:
        child_depth = get_max_deepth(child)
        max_child_depth = max(max_child_depth, child_depth)
    
    return max_child_depth + 1


def compute_tip_to_root_branch_length_variance(tree: object) -> int:
    """
    Compute the variance of branch lengths from tips to root.

    Args:
        tree (object): Tree object.

    Returns:
        int: Variance of branch lengths.
    """
    tip_to_root_branch_lengths = []
    for leaf in tree.iter_leaves():
        branch_length = tree.get_distance(leaf)
        tip_to_root_branch_lengths.append(branch_length)
    branch_length_variance = 0
    if len(tip_to_root_branch_lengths) > 1:
        variance=float(np.var(tip_to_root_branch_lengths))
        branch_length_variance = variance
    return branch_length_variance

def calculate_species_num(node: object) -> int:
    """
    Calculate the number of unique species under a node.

    Args:
        node (object): Tree node.

    Returns:
        int: Number of unique species.
    """
    species_num=len(get_species_set(node))
    return species_num
    
def sps_dup_num(sps_list:list, unique_sps:list)->int:
    sps_num_dic = {i: 0 for i in unique_sps}
    sps_dups = set()

    for sps in sps_list:
        if sps in sps_num_dic:
            sps_num_dic[sps] += 1
            if sps_num_dic[sps] > 1:
                sps_dups.add(sps)

    return len(sps_dups)

def calculate_gd_num(Phylo_t: object) -> int:
    """
    Calculate the number of gene duplication events in a phylogenetic tree.

    Args:
        Phylo_t (object): Phylogenetic tree object.

    Returns:
        int: Number of gene duplication events.
    """
    gd_num=0
    gd_node_names=find_dup_node_simple(Phylo_t)
    for node in gd_node_names:
        clade=node
        sps=[leaf.split('_')[0] for leaf in clade.get_leaf_names()]
        unique_sps=set(sps)
        if len(unique_sps) >5:
            if sps_dup_num(sps,unique_sps) > len(unique_sps)*0.2:
                gd_num+=1
        else:
            if sps_dup_num(sps,unique_sps) >=1:
                gd_num+=1
    
    return gd_num

def count_common_elements(scubclade1_sps: set, subclade2_sps: set) -> int:
    """
    Count the number of species that appear more than once in the list.

    Args:
        scubclade1_sps (list): List of subclade1 species names.
        scubclade2_sps (list): List of subclade2 species names.

    Returns:
        int: Number of duplicated species.
    """
   
    return len(scubclade1_sps & subclade2_sps)

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
    
def find_tre_dup(Phylo_t: object) -> list:
    """
    Find all duplication events in a phylogenetic tree and return their sequence relationships.

    Args:
        Phylo_t (object): Phylogenetic tree object.

    Returns:
        tuple: A list of duplication event relationships (as strings) and a set of gene family leaf names.
    """
    tre_ParaL = []
    if not hasattr(Phylo_t, "get_leaf_names") or not hasattr(Phylo_t, "get_descendant_evol_events"):
        raise ValueError("Input object does not have required tree methods.")
    GF_leaves_S = set(Phylo_t.get_leaf_names())
    events = Phylo_t.get_descendant_evol_events()
    for ev in events:
        if ev.etype == "D":
            tre_ParaL.append(",".join(ev.in_seqs) + "<=>" + ",".join(ev.out_seqs))
    return tre_ParaL, GF_leaves_S


def realign_branch_length(Phylo_t1:object)->object:
    Phylo_t1.ladderize()
    Phylo_t1.resolve_polytomy(recursive=True)
    Phylo_t1.sort_descendants("support")
    max_deep=get_max_deepth(Phylo_t1)
    for node in Phylo_t1.traverse():
        if not node.is_root():
            node.dist=1
            degree=node.get_distance(node.get_tree_root()) + 1
            deep=get_max_deepth(node)
            node.dist=max_deep-deep-degree
    clade_up=Phylo_t1.get_children()[0]
    clade_down=Phylo_t1.get_children()[1]
    difference=abs(get_max_deepth(clade_up)-get_max_deepth(clade_down))+1
    clade_up.dist=clade_up.dist+difference  
    clade_down.dist=clade_down.dist+difference   
    
    return Phylo_t1

def rejust_root_dist(sptree):
    clade_up=sptree.get_children()[0]
    clade_down=sptree.get_children()[1]
    if len(clade_up)>len(clade_down):
        clade_up.dist=1 
        if clade_down.is_leaf():
            clade_down.dist=get_max_deepth(sptree)-1
        else:
            clade_down.dist=get_max_deepth(sptree)-get_max_deepth(clade_down)
    else:
        clade_down.dist=1 
        if clade_up.is_leaf():
            clade_up.dist=get_max_deepth(sptree)-1
        else:
            clade_up.dist=get_max_deepth(sptree)-get_max_deepth(clade_up)

    return sptree

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

def find_dup_node(
    gene_tree: object,
    species_tree: object,
    gd_support: int = 50,
    clade_support: int = 0,
    dup_species_num: int = 2,
    dup_species_percent: float = 0.2,
    max_topology_distance: int = 1
) -> list:
    """
    Find duplication nodes using reconciliation events,
    but return them in strict post-order (tree2gd-compatible).
    """

    # 1️⃣ 先从 event 中拿到 duplication 对应的节点（集合）
    dup_nodes = set()

    events = gene_tree.get_descendant_evol_events()
    for event in events:
        if event.etype != "D":
            continue

        node_names = list(event.in_seqs) + list(event.out_seqs)
        ca = gene_tree.get_common_ancestor(node_names)

        dup_nodes.add(ca)

    # 2️⃣ 再 post-order 遍历，按顺序筛选
    dup_node_list = []

    for node in gene_tree.traverse("postorder"):
        if node not in dup_nodes:
            continue

        if not judge_support(node.support, gd_support):
            continue

        children = node.get_children()
        if len(children) != 2:
            continue
        child_a, child_b = children

        if child_a.support < clade_support or child_b.support < clade_support:
            continue

        species_set = get_species_set(node)

        mapped = map_gene_tree_to_species(species_set, species_tree)
        node.add_feature('map', mapped.name)

        if len(species_set) == 1:
            dup_node_list.append(node)
            continue

        dup_sps = count_common_elements(
            get_species_set(child_a),
            get_species_set(child_b)
        )
        dup_percent = dup_sps / len(species_set)

        if dup_sps < dup_species_num or dup_percent < dup_species_percent:
            continue

        mapped_a = map_gene_tree_to_species(get_species_set(child_a), species_tree)
        mapped_b = map_gene_tree_to_species(get_species_set(child_b), species_tree)

        if species_tree.get_distance(
            mapped_a, mapped_b, topology_only=True
        ) > max_topology_distance:
            continue

        dup_node_list.append(node)

    return dup_node_list


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
