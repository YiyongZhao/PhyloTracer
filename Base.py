import pandas as pd
from ete3 import PhyloTree
import random
import numpy as np

def generate_sps_voucher(sps_num:int) -> list:
    characters = [chr(i) for i in range(65, 91)] + [chr(i) for i in range(97, 123)] + [str(i) for i in range(10)]
    unique_strings = set()

    while len(unique_strings) < sps_num:
        unique_strings.add(''.join(random.sample(characters, 3)))

    return sorted(list(unique_strings))


def gene_id_transfer(gene2taxa_list:str) -> dict:
    gene2taxa_dic=read_and_return_dict(gene2taxa_list)
    taxa_list=list(set(gene2taxa_dic.values()))
    taxa2voucher_dic=dict(zip(taxa_list,generate_sps_voucher(len(taxa_list))))
    voucher2taxa_dic=dict(zip(taxa2voucher_dic.values(),taxa2voucher_dic.keys()))
    gene2new_named_gene_dic = {gene: f"{taxa2voucher_dic[gene2taxa_dic[gene]]}_{gene}" for gene in gene2taxa_dic.keys()}
    new_named_gene2gene_dic=dict(zip(gene2new_named_gene_dic.values(),gene2new_named_gene_dic.keys()))
    return gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic
#gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("gene2taxa.list")

def read_and_return_dict(filename, separator="\t") -> dict:
    df=pd.read_csv(filename,sep=separator,header=None)
    return df.set_index([0])[1].to_dict()

def rename_input_tre(Phylo_t:object, gene2new_named_gene_dic:dict) -> object:
    for node in Phylo_t.traverse():
        if node.name in gene2new_named_gene_dic.keys():
            node.name = gene2new_named_gene_dic[node.name]
    return Phylo_t

def read_tree(tre_path:str) -> object:
    return PhyloTree(tre_path)

######################################################################################################################
def isRoot_Judger(Phylo_t:object)->bool:#Translate the function that determines whether the input phylogenetic tree has a negated root.
    if len(Phylo_t.get_children()) ==2:
        return True

def root_tre_with_midpoint_outgroup(Phylo_t:object)->object:#Rooting the phylogenetic tree using the midpoint outgroup method.
    mid_node=Phylo_t.get_midpoint_outgroup()
    Phylo_t.set_outgroup(mid_node)
    return Phylo_t

def num_tre_node(Phylo_t:object)->object:#Numbering the nodes in the tree.
    i = 1
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            node.name = "node" + str(i)
            i += 1
    return Phylo_t

def get_species_list(Phylo_t:object)->list:
    leaves = Phylo_t.get_leaf_names()
    species_lst = []
    for leaf in leaves:
        species = leaf.split("_")[0]  # Assuming that the species name is located in the first part of the node name (leaf) with “_” as the delimiter
        species_lst.append(species)
    return species_lst

class TreeNode:
    def __init__(self, val, branch_length):
        self.val = val
        self.branch_length = branch_length
        self.children = []

def get_max_deepth(root:object)->int:
    if not root:
        return 0
    
    max_child_depth = 0
    for child in root.children:
        child_depth = get_max_deepth(child)
        max_child_depth = max(max_child_depth, child_depth)
    
    return max_child_depth + 1

def compute_tip_to_root_branch_length_variance(tree:object)->int:
    tip_to_root_branch_lengths = []
    for leaf in tree.iter_leaves():
        branch_length = tree.get_distance(leaf)
        tip_to_root_branch_lengths.append(branch_length)

    branch_length_variance = 0
    if len(tip_to_root_branch_lengths) > 1:
        variance=float(np.var(tip_to_root_branch_lengths))
        branch_length_variance = variance

    return branch_length_variance

def calculate_species_overlap(gene_tree:object)->int:
    up_clade=gene_tree.children[1]
    down_clade=gene_tree.children[0]
    species_list_a = get_species_list(up_clade)
    species_list_b = get_species_list(down_clade)
    overlap_ratio = len(set(species_list_a) & set(species_list_b)) / len(set(species_list_a) | set(species_list_b))

    return overlap_ratio

def calculate_species_num(node:object)->int:# Obtain the number of species under a node
    leaf_list=node.get_leaf_names()
    species_num=len(set(i.split('_')[0] for i in leaf_list))    
    return species_num

