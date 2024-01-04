from Ortho_Retriever import *
from __init__ import *

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

######################################################################################################################
def get_RF_dic(voucher2taxa_dic,trees:list,gene2new_named_gene_dic:dict,tre_ID:str,tre_path:str,renamed_len_dic:dict,new_named_gene2gene_dic:dict,sptree:object)->dict:
    sptree=rename_species_tree(sptree, voucher2taxa_dic)
    RF_dic={}
    
    def calculate_RF_distance(Phylo_t_OG_L:object,sptree:object)->int:# RF comparison between each individual gene tree and the species tree
        for leaf in Phylo_t_OG_L:
            leaf.name=leaf.name[0:3]# Processing leaf names, such as changing ABC_001 to ABC, in order to compare with the species tree
        RF=Phylo_t_OG_L.robinson_foulds(sptree)[0]
        return RF

    processed_trees = set()
    
    for i in trees:
        if i in processed_trees:
            continue
        #if not yes_or_no_root(i):
             #mid_node=i.get_midpoint_outgroup()
            #i.set_outgroup(mid_node)
        principal_gene_S,filtered_offcut_ev_seqs_L0=offcut_tre(i,renamed_len_dic)
        minor_othologs_L=[]
        minor_othologs_L=iterator(filtered_offcut_ev_seqs_L0,i,gene2new_named_gene_dic,minor_othologs_L,tre_path,renamed_len_dic)
        ordered_name_OG_L=rename_OGs_tre_name(principal_gene_S,minor_othologs_L,tre_ID)
        RF=0
        for tre_name,OG_S in ordered_name_OG_L: 
            OG_L = [new_named_gene2gene_dic[OG] for OG in OG_S]
            Phylo_t0=read_phylo_tree(tre_path)
            Phylo_t=root_tre_with_midpoint_outgroup(Phylo_t0)
            Phylo_t_OG_L=extract_tree(OG_L,Phylo_t)
            Phylo_t_OG_L=rename_input_tre(Phylo_t_OG_L,gene2new_named_gene_dic)
            #o.write(tre_name+"\t"+Phylo_t_OG_L.write()+"\n")
            RF+=calculate_RF_distance(Phylo_t_OG_L,sptree)
        RF_dic[i]=RF
        processed_trees.add(i)
    return RF_dic
######################################################################################################################
def is_same_species(node:object)->bool:# Determining if the genes under a node are from the same species
    if calculate_species_num(node) >1:
        return True
                
def traverse(node:object)->str: #Iterate through the upward and downward branches of a node
    if is_same_species(node):
        return node.name+'-'+traverse(node.children[0])+'-'+traverse(node.children[1])
    return node.name

def get_reroot_list(Phylo_t0:object,node_name_list:list)->list:#Get the list of new trees obtained by re-rooting a multi-copy tree with each node as an outgroup
    re_root_list=[]
    for i in node_name_list:
        Phylo_t1=Phylo_t0.copy('newick')
        node=Phylo_t1&i
        Phylo_t1.set_outgroup(node)
        Phylo_t1.name=i
        re_root_list.append(Phylo_t1)
    return re_root_list

def manager(Phylo_t0:object)->list:#Get a list of tree objects after rooting them
    Phylo_t0=num_tre_node(Phylo_t0)
    node_name_list=traverse(Phylo_t0).split('-')#Get the node_name_list of the tree that needs to be re-rooted
    node_name_list.remove('N1')
    root_list=get_reroot_list(Phylo_t0,node_name_list)#The list of tree objects that require root reassignment needs to be translated as follows:
    return root_list

def rename_output_tre(Phylo_t:object, new_named_gene2gene_dic:dict,tre_ID:str) -> object:#Translate the renamed phylogenetic tree back to the original tree.
    for node in Phylo_t.traverse():
        if node.name in new_named_gene2gene_dic.keys():
            node.name = new_named_gene2gene_dic[node.name]
    return Phylo_t.write(outfile=str(tre_ID)+'.nwk')

def root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic,sptree,voucher2taxa_dic):
    for tre_ID,tre_path in tre_dic.items():
        Phylo_t0=read_phylo_tree(tre_path)
        Phylo_t1=root_tre_with_midpoint_outgroup(Phylo_t0)
        if len(set(get_species_list(Phylo_t1))) <= 2:
            Phylo_t1.write(outfile=str(tre_ID)+'.nwk')
            continue
        Phylo_t2=rename_input_tre(Phylo_t1,gene2new_named_gene_dic)
        root_list=manager(Phylo_t2)
        tre_ParaL,GF_leaves_S = find_tre_dup(Phylo_t2)
        if len(tre_ParaL) ==None:#Find gene duplications, calculate the excess depth of outgroups without duplications, and select the tree with the maximum excess depth of outgroups.
            deeped_trees= filter_max_deep(root_list)
            if len(deeped_trees) > 1:
                vared_trees=filter_min_var(deeped_trees)
                if len(vared_trees) >1:
                    RF_dic=get_RF_dic(voucher2taxa_dic,vared_trees,gene2new_named_gene_dic,tre_ID,tre_path,renamed_len_dic,new_named_gene2gene_dic,sptree)   
                    RFed_trees=filter_min_RF(RF_dic)
                    rename_output_tre(RFed_trees[0],new_named_gene2gene_dic,tre_ID)
                else:
                    rename_output_tre(vared_trees[0],new_named_gene2gene_dic,tre_ID) 
            else:
                rename_output_tre(deeped_trees[0],new_named_gene2gene_dic,tre_ID)
        else:
            GDed_trees=filter_min_GD(root_list)#Identify gene duplications and if duplications exist, calculate the number of duplications to filter, and select the tree with the minimum number of duplications
            if len(GDed_trees)>1:
                max_overlaped_trees=filter_max_overlap(GDed_trees)
                if len(max_overlaped_trees) >1:
                        vared_trees=filter_min_var(max_overlaped_trees)
                        if len (vared_trees) > 1:
                            RF_dic=get_RF_dic(voucher2taxa_dic,vared_trees,gene2new_named_gene_dic,tre_ID,tre_path,renamed_len_dic,new_named_gene2gene_dic,sptree)   
                            RFed_trees=filter_min_RF(RF_dic) 
                            rename_output_tre(RFed_trees[0],new_named_gene2gene_dic,tre_ID) 
                        else:
                            rename_output_tre(vared_trees[0],new_named_gene2gene_dic,tre_ID) 
                else:
                    rename_output_tre(max_overlaped_trees[0],new_named_gene2gene_dic,tre_ID)
            else:
                rename_output_tre(GDed_trees[0],new_named_gene2gene_dic,tre_ID)

if __name__ == "__main__":
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    len_dic=read_and_return_dict('length_list')
    renamed_len_dic=rename_len_dic(len_dic,gene2new_named_gene_dic)
    sptree=PhyloTree('sptree')
    tre_dic=read_and_return_dict('GF_list')   
    root_main(tre_dic, gene2new_named_gene_dic, renamed_len_dic, new_named_gene2gene_dic,sptree,voucher2taxa_dic)
