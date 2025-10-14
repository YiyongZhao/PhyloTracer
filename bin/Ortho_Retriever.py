from __init__ import *

def rename_len_dic(len_dic:dict, gene2new_named_gene_dic:dict) -> dict:
    """
    Rename the keys of the length dictionary using the gene2new_named_gene_dic mapping.
    Args:
        len_dic (dict): Dictionary mapping gene names to lengths.
        gene2new_named_gene_dic (dict): Mapping from original gene names to new names.
    Returns:
        dict: Dictionary with renamed keys and original values.
    """
    return {gene2new_named_gene_dic[key]: value for key,value in len_dic.items() if key in gene2new_named_gene_dic.keys()}

def count_sps_num(ev_seqs:set)->int:
    """
    Count the number of unique species in a set of gene names.
    Args:
        ev_seqs (set): Set of gene names.
    Returns:
        int: Number of unique species (by prefix).
    """
    sps_set={gene[0:3] for gene in ev_seqs }
    return len(sps_set) 
  
def extract_tree(ev_seqs:set,Phylo_t:object)-> object:
    """
    Prune the phylogenetic tree to keep only the specified gene set.
    Args:
        ev_seqs (set): Set of gene names to retain.
        Phylo_t (object): Phylogenetic tree object.
    Returns:
        object: Pruned phylogenetic tree.
    """
    Phylo_t.prune(ev_seqs)
    return Phylo_t

def offcut_tre(Phylo_t:object,renamed_len_dic:dict)-> list:
    """
    Identify and separate offcut gene sets from the phylogenetic tree based on length and species count.
    Args:
        Phylo_t (object): Phylogenetic tree object.
        renamed_len_dic (dict): Mapping from gene names to lengths (renamed).
    Returns:
        tuple: (principal_gene_S, filtered_offcut_ev_seqs_L0)
            principal_gene_S (set): Main gene set after offcut removal.
            filtered_offcut_ev_seqs_L0 (list): Filtered list of offcut gene sets.
    """
    tre_ParaL,GF_leaves_S = find_tre_dup(Phylo_t)
    offcut_ev_seqs_L0=[]
    for ev_seqs in tre_ParaL:
        fd=ev_seqs.strip().split("<=>")
        ev_seqs1,ev_seqs2=set(fd[0].split(",")),set(fd[1].split(","))
        sps_num1,seq_num1,sps_num2,seq_num2 = count_sps_num(ev_seqs1),len(ev_seqs1),count_sps_num(ev_seqs2),len(ev_seqs2)
        if sps_num1 == sps_num2:
            if seq_num1 > seq_num2:
                offcut_ev_seqs_L0.append(ev_seqs2)
            elif seq_num1 < seq_num2:
                offcut_ev_seqs_L0.append(ev_seqs1)
            elif seq_num1 == seq_num2:
                ev_seqs1_avg_len,ev_seqs2_avg_len = sum([renamed_len_dic[Gene] for Gene in ev_seqs1 ])/sps_num1,sum([renamed_len_dic[Gene] for Gene in ev_seqs2 ])/sps_num2
                if ev_seqs1_avg_len >= ev_seqs2_avg_len:
                    offcut_ev_seqs_L0.append(ev_seqs2)
                else:
                    offcut_ev_seqs_L0.append(ev_seqs1)
        else:
            if sps_num1 > sps_num2:
                offcut_ev_seqs_L0.append(ev_seqs2)
            else:
                offcut_ev_seqs_L0.append(ev_seqs1)
    offcut_Gene_S=set()
    for ev_seqs in offcut_ev_seqs_L0:
        offcut_Gene_S = offcut_Gene_S | ev_seqs
    principal_gene_S=GF_leaves_S - offcut_Gene_S
    filtered_offcut_ev_seqs_L0=[]
    for ev_seqs in offcut_ev_seqs_L0:
        sps_num = count_sps_num(ev_seqs)
        if sps_num >= 2:
            filtered_offcut_ev_seqs_L0.append(ev_seqs)
    filtered_offcut_ev_seqs_L0=rm_dup(filtered_offcut_ev_seqs_L0)
    return principal_gene_S,filtered_offcut_ev_seqs_L0# filtered_offcut_ev_seqs_L0=offcut_ev_seqs_L0

def rm_dup(paralogs_L): #[{1,2},{2,3}]
    for ev_seqs1 in paralogs_L:
        for ev_seqs2 in paralogs_L:
            if ev_seqs1.issubset(ev_seqs2):
                if ev_seqs1 not in paralogs_L:
                    paralogs_L.remove(ev_seqs1)
            elif ev_seqs2.issubset(ev_seqs1):
                if ev_seqs2 not in paralogs_L:
                    paralogs_L.remove(ev_seqs2)
    return paralogs_L

def split_offcut_ev_seqs(offcut_ev_seqs_L0:list) ->list:  #[{1,2},{2,3}]
    othologs_L=[]
    paralogs_L=[]
    for ev_seqs in offcut_ev_seqs_L0:
        sps_num,seq_num = count_sps_num(ev_seqs),len(ev_seqs)
        if sps_num >=2:
            if seq_num > sps_num:
                paralogs_L.append(ev_seqs)
            else:
                othologs_L.append(ev_seqs)
    paralogs_L=rm_dup(paralogs_L)
    return othologs_L,paralogs_L

def iterator(offcut_ev_seqs_L0:list,Phylo_t:object,new_named_gene2gene_dic:dict,minor_othologs_L:list,tre_path:str,renamed_len_dic:dict) ->list:
    othologs_L,paralogs_L=split_offcut_ev_seqs(offcut_ev_seqs_L0)
    minor_othologs_L+=(othologs_L)
    if paralogs_L !=[]:    
        for i,ev_seqs in enumerate(paralogs_L):
            Phylo_t=read_phylo_tree(tre_path)
            if is_rooted(Phylo_t):
                pass
            else:
                Phylo_t=root_tre_with_midpoint_outgroup(Phylo_t)
            Phylo_t=rename_input_tre(Phylo_t,new_named_gene2gene_dic)
            Phylo_t = extract_tree(ev_seqs,Phylo_t)
            principal_gene_S,offcut_ev_seqs_L0=offcut_tre(Phylo_t,renamed_len_dic)
            minor_othologs_L.append(principal_gene_S)
            iterator(offcut_ev_seqs_L0,Phylo_t,new_named_gene2gene_dic,minor_othologs_L,tre_path,renamed_len_dic)
    return minor_othologs_L
      
def rename_OGs_tre_name(principal_gene_S:list,minor_othologs_L:list,tre_ID:str)->list:
    """
    Generate ordered names for ortholog groups and associate them with gene sets.
    Args:
        principal_gene_S (list): Main gene set.
        minor_othologs_L (list): List of minor ortholog sets.
        tre_ID (str): Tree identifier.
    Returns:
        list: List of tuples (tree_name, gene_set).
    """
    OG_count = len(minor_othologs_L)+1
    tre_name_L =["T"+str(i+1) for i in range(OG_count)]
    sps_num_L=[len(OG) for OG in minor_othologs_L]
        
    tre_name2sps_num_dict=dict(zip(tre_name_L,sps_num_L))
    tre_name2OG_dict=dict(zip(tre_name_L,minor_othologs_L))
    orderd_OG_sps_num_L = sorted(tre_name2sps_num_dict.items(), key=lambda x: x[1],reverse =True)
        
    ordered_name_OG_L=[]
    ordered_name_OG_L.append((str(tre_ID)+"_"+tre_name_L[0]+"_"+str(len(principal_gene_S)),principal_gene_S))
    index=0
    for name,count in orderd_OG_sps_num_L:
        index+=1
        ordered_name_OG_L.append((str(tre_ID)+"_"+tre_name_L[index]+"_"+str(count),tre_name2OG_dict[name]))
    return ordered_name_OG_L

def get_single_copy_trees(Phylo_t1:object,renamed_len_dic:dict,gene2new_named_gene_dic:dict,new_named_gene2gene_dic:dict,tre_path:str,tre_ID:str)->list:
    """
    Extract single-copy ortholog trees from the phylogenetic tree.
    Args:
        Phylo_t1 (object): Phylogenetic tree object.
        renamed_len_dic (dict): Mapping from gene names to lengths (renamed).
        gene2new_named_gene_dic (dict): Mapping from gene names to new names.
        new_named_gene2gene_dic (dict): Mapping from new gene names to original gene names.
        tre_path (str): Path to the tree file.
        tre_ID (str): Tree identifier.
    Returns:
        list: List of tuples (tree_name, pruned_tree_object).
    """
    trees=[]
    principal_gene_S,filtered_offcut_ev_seqs_L0=offcut_tre(Phylo_t1,renamed_len_dic)
    minor_othologs_L=[]
    minor_othologs_L=iterator(filtered_offcut_ev_seqs_L0,Phylo_t1,gene2new_named_gene_dic,minor_othologs_L,tre_path,renamed_len_dic)
    ordered_name_OG_L=rename_OGs_tre_name(principal_gene_S,minor_othologs_L,tre_ID)
    for tre_name,OG_S in ordered_name_OG_L:

        OG_L = [new_named_gene2gene_dic[OG] for OG in OG_S]
        Phylo_t0=read_phylo_tree(tre_path)
        if is_rooted(Phylo_t0):
            Phylo_t=Phylo_t0
        else:
            Phylo_t=root_tre_with_midpoint_outgroup(Phylo_t0)
        Phylo_t_OG_L=extract_tree(OG_L,Phylo_t)
        trees.append((tre_name,Phylo_t_OG_L))
    return trees

#ordered_name_OG_L=rename_OGs_tre_name(principal_gene_S,minor_othologs_L,tre_ID)
##########################################################################
def split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic):
    """
    Main function to split gene family trees into single-copy ortholog trees and write results.
    Args:
        tre_dic (dict): Mapping from tree IDs to tree file paths.
        gene2new_named_gene_dic (dict): Mapping from gene names to new names.
        new_named_gene2gene_dic (dict): Mapping from new gene names to original gene names.
        renamed_len_dic (dict): Mapping from gene names to lengths (renamed).
    Returns:
        None
    """
    o=open('ortho_retriever_summary.txt','w')
    o.write('tre_name'+'\t'+'single_copy_tree'+'\n')
    processed_lines=[]
    for tre_ID,tre_path in tre_dic.items():
        Phylo_t0=read_phylo_tree(tre_path)
        if is_rooted(Phylo_t0):
            Phylo_t1=Phylo_t0
        else:
            Phylo_t1=root_tre_with_midpoint_outgroup(Phylo_t0)

        Phylo_t2=rename_input_tre(Phylo_t1,gene2new_named_gene_dic)
        trees=get_single_copy_trees(Phylo_t2,renamed_len_dic,gene2new_named_gene_dic,new_named_gene2gene_dic,tre_path,tre_ID)
        for clade in trees:
            tre_name=clade[0]
            Phylo_t_OG_L=clade[1]
            processed_lines.append(tre_name + "\t" + Phylo_t_OG_L.write()+ "\n")
    sorted_lines = sorted(processed_lines, key=lambda x: len(x.split('\t')[1]), reverse=True)
    for line in sorted_lines:
        o.write(line)
    o.close()

if __name__ == "__main__":
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    len_dic=read_and_return_dict('length')
    renamed_len_dic=rename_len_dic(len_dic,gene2new_named_gene_dic)
    tre_dic=read_and_return_dict('GF_list.txt')   
    split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic,renamed_len_dic)
