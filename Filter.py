from ete3 import PhyloTree
import numpy as np
from Base import *
from Ortho_Split import *
from Phylo_Rooting import *

def filter_min_RF(RF_dic:dict)->list:#Filter RF and select the tree with the minimum value
    RF_trees=RF_dic.keys()
    RF_list=RF_dic.values()
    min_RF=min(RF_list)
    RFed_trees = [RF_tree for RF_tree, RF in zip(RF_trees, RF_list) if RF == min_RF]
    
    return RFed_trees

def filter_max_deep(RFed_trees:list)->list:#Filter outgroups and take the maximum depth
    deep_list=[]
    for i in RFed_trees:
        up_clade=i.children[1]
        down_clade=i.children[0]
        if len(up_clade.get_leaf_names())>len(down_clade.get_leaf_names()):
            
        #if not h.is_leaf():
            t=get_max_deepth(down_clade)
        else:
            t=get_max_deepth(up_clade)
        deep_list.append(t)
    max_depth = max(deep_list)
    deeped_trees = [RFed_tree for RFed_tree, depth in zip(RFed_trees, deep_list) if depth == max_depth]
    
    return deeped_trees

def filter_min_var(deeped_trees:list)->list:#Filter branches by variance and select the tree with the minimum value
    var_list=[]
    for i in deeped_trees:
        up_clade=i.children[1]
        down_clade=i.children[0]
        if len(up_clade.get_leaf_names())>len(down_clade.get_leaf_names()):
            
        #if not h.is_leaf():
            var=compute_tip_to_root_branch_length_variance(up_clade)
        else:
            var=compute_tip_to_root_branch_length_variance(down_clade)
        var_list.append(var)
    min_var=min(var_list)
    vared_trees = [deeped_tree for deeped_tree, var in zip(deeped_trees, var_list) if var == min_var]
    
    return vared_trees
   
def filter_min_GD(RFed_trees:list)->list:#Filter by GD_num and select the tree with the minimum value
    GD_list=[]
    for i in RFed_trees:
        tre_ParaL,GF_leaves_S = find_tre_dup(i) 
        GD_num=len(tre_ParaL)
        GD_list.append(GD_num)
    min_GD=min(GD_list)
    GDed_trees = [RFed_tree for RFed_tree, GD_num in zip(RFed_trees, GD_list) if GD_num == min_GD]
    return  GDed_trees  

def filter_max_overlap(GDed_trees:list)->list:#Filter by species coverage and select the tree with the maximum value
    overlap_list=[]
    for i in GDed_trees:
        overlap_ratio=calculate_species_overlap(i)
        overlap_list.append(overlap_ratio)
    max_overlap_ratio=max(overlap_list)
    max_overlaped_trees = [GDed_tree for GDed_tree, overlap_ratio in zip(GDed_trees, overlap_list) if overlap_ratio == max_overlap_ratio]
    return max_overlaped_trees 
