from __init__ import *
import re



def count_common_elements(set1, set2):
   
    return len(set1.intersection(set2))

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
            mapped_a = map_gene_tree_to_species(get_species_set(child_a), species_tree)
            mapped_b = map_gene_tree_to_species(get_species_set(child_b), species_tree)
            dup_sps=count_common_elements(get_species_set(child_a), get_species_set(child_b))
            dup_percent = dup_sps / len(get_species_set(common_ancestor_node))
            if common_ancestor_node.support>=gd_support:
                if child_a.support >= clade_support and child_b.support >= clade_support:
                    if dup_percent>=dup_species_percent:
                        if species_tree.get_distance(mapped_a, mapped_b, topology_only=True) <= max_topology_distance:
                            dup_node_list.append(common_ancestor_node)

                       
    return dup_node_list
    
def get_maptree_internal_node_name_set(node,sptree):
    sps=get_species_set(node)

    node1=node.up
    map_nodename1=get_maptree_name(sptree,get_species_set(node1))
    map_node1=sptree&map_nodename1

    names=[]
    for i in sps:
        clade=sptree&i
        path=get_two_nodes_path_str(clade,map_node1)
        names+=path

    return set(names)


def get_two_subclade_maptree_node_name_lst(max_clade,sptree):
    clade_up=max_clade.get_children()[0]
    clade_up_set=get_maptree_internal_node_name_set(clade_up,sptree)
    clade_down=max_clade.get_children()[1]
    clade_down_set=get_maptree_internal_node_name_set(clade_down,sptree)
    up_down_lst=list(clade_up_set)+list(clade_down_set)
    return up_down_lst 

def get_maptree_internal_node_name_count_dic(max_clade,max_clade2sp,sptree):
    up_down_lst=get_two_subclade_maptree_node_name_lst(max_clade,sptree)
    dic={i.name:0 for i in max_clade2sp.traverse()}
    for i in up_down_lst :
        if i in dic :
            dic[i]+=1
    keys_with_zero_value = [key for key, value in dic.items() if value == 0 and key not in sptree.get_leaf_names()]
    return dic,keys_with_zero_value

def get_maptree_name(sptree, sps_set):
    if len(sps_set)!=1:
        com = sptree.get_common_ancestor(sps_set)
        return com.name
    else:
        com=sptree&list(sps_set)[0]
        return com.name
    
def get_two_nodes_path_str(start_node, end_node):
    nodes = []
    current_node = start_node
    while current_node != end_node:
        nodes.append(current_node.name)
        current_node = current_node.up
        if current_node is None:
            break
    nodes.append(end_node.name)
    re_nodes=list(reversed(nodes))
   
    return re_nodes

def get_tips_to_clade_path_lst(taget_node:object,dic)->list:#taget_node也就是gd node
    path_str_lst=[]
    for i in taget_node:
        path_str=get_two_nodes_path_str(i,taget_node)
        numlist=[dic[k] for k in path_str]
        result = '->'.join([f"{name}({count})" for name, count in zip(path_str, numlist)])

        path_str_lst.append(result)
    return path_str_lst
        
def get_maptree_node_count_dic(sp_list,map_clade):#{'S1':2,'A':2}
    sp_dic={i.name:0 for i in map_clade.traverse()}
    for i in sp_list:
        if i in sp_dic:
            sp_dic[i]+=1
    
    tips=set([j for j in sp_dic.keys() if not j.startswith('S') ])
    for k,v in sp_dic.items():
        if k.startswith('S'):
            if v==0:
                t=map_clade&k
                sp=get_species_set(t)
                if len(sp&tips)!=0:
                    sp_dic[k]=1
            elif v==1:
                t=map_clade&k
                sp=get_species_set(t)
                if len(sp&tips)!=0:
                    jiao=list(sp&tips)[0]
                    if sp_dic[jiao]==2:
                        sp_dic[k]=2   
    return sp_dic
    
            
def add_extra_maptree_node_name(lst,sptree):
    lst1=lst.copy()
    for i in sptree.traverse():
        if not i.is_leaf():
            if not i.name in lst1:
                if len(set(lst)&get_species_set(i))!=0:
                    lst1.add(i.name)       
    return lst1

def get_gene_loss_list(clade_up,clade_down,dup_sps):
    gene_loss=[]
    for leaf in clade_up.get_leaf_names():
        if leaf.split('_')[0]  in dup_sps:
            continue
        else: 
            gene_loss.append(leaf)

    for leaf2 in clade_down.get_leaf_names():
        if leaf2.split('_')[0] in dup_sps:
            continue
        else: 
            gene_loss.append(leaf2)

    return gene_loss

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

def get_path_str_with_count_num_lst(tre_id,gd_id,genetree,renamed_sptree,out,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic,path2_treeid_dic):#统计一颗树中所有gd下不同的丢失情况
    path_str_with_count_num_lst = []
    dup_node_name_list = find_dup_node(genetree,renamed_sptree)
    for i in dup_node_name_list:
        sp = get_species_set(i)
        max_clade2sp = mapp_gene_tree_to_species(sp, renamed_sptree)
        clade_up = i.get_children()[0]
        clade_down = i.get_children()[1]
        new_path_str_lst = []
        if len(sp) != 1:
            dic,keys_with_zero_value=get_maptree_internal_node_name_count_dic(i,max_clade2sp,renamed_sptree)
            path_str_lst = get_tips_to_clade_path_lst(max_clade2sp,dic)
            for i1 in path_str_lst:
                if i1 in path2_treeid_dic:
                    path2_treeid_dic[i1].append(f'{tre_id}-{gd_id}')
                else:
                    path2_treeid_dic[i1]=[f'{tre_id}-{gd_id}']

            path_str_with_count_num_lst+=path_str_lst

           
            for path_str in path_str_lst:
                parts = path_str.split('->')
                last = parts[-1]
                key = last.split('(')[0]
                if key in voucher2taxa_dic:
                    parts[-1] = voucher2taxa_dic[key] + last[len(key):]
                new_path_str_lst.append('->'.join(parts))
        
        species_set=get_species_list(max_clade2sp)
        appeared_species = set()
        gene_pairs = gene_pair(i)
        for gene_pair_str in gene_pairs:
            out.write(str(tre_id) + '\t' +str(gd_id) + '\t'+ str(i.support)+'\t'+voucher2taxa_dic.get(max_clade2sp.name,max_clade2sp.name) +'\t')
            gene_a, gene_b = gene_pair_str.split('-')
            if gene_a.split('_')[0] == gene_b.split('_')[0]:
                species_code = gene_a.split('_')[0]
            else:
                if gene_a == 'null':
                    species_code = gene_b.split('_')[0]
                else:
                    species_code = gene_a.split('_')[0]
            appeared_species.add(species_code)
            out.write(str(voucher2taxa_dic[species_code]) + '\t' +
                        str(new_named_gene2gene_dic.get(gene_a, 'null')) + '\t' +
                        str(new_named_gene2gene_dic.get(gene_b, 'null')) + '\t' +
                        '@'.join(new_path_str_lst)+'\n'
            )

        for s in species_set:
            if s not in appeared_species:
                out.write(str(tre_id) + '\t' + str(gd_id) + '\t'+ str(i.support)+'\t' + voucher2taxa_dic.get(s,s)  + '\tNA\tNA\tNA\t'+ '@'.join(new_path_str_lst) + '\n')
        gd_id+=1
        
    set_path_str_with_count_num_lst = path_str_with_count_num_lst.copy()

    
    result = sorted(set_path_str_with_count_num_lst, reverse=True)
    
    return result,gd_id

def num_sptree(sptree):
    n=0
    for i in sptree.traverse('postorder'):
        if not i.is_leaf():
            i.name='S'+str(n)
            n+=1
    sptree.sort_descendants()
    sptree.write(outfile='numed_sptree.nwk',format=1)
    return sptree

def split_dict_by_first_last_char(original_dict):
    split_dicts = {}

    for key, value in original_dict.items():
        first_char = key.split('->')[0].split('(')[0]
        last_char = key.split('->')[-1].split('(')[0]


        new_key = f"{first_char}_{last_char}"

        if new_key not in split_dicts:
            split_dicts[new_key] = {}

        split_dicts[new_key][key] = value

    return split_dicts

def get_path_str_num_dic(tre_dic,sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic):
    renamed_sptree=rename_input_tre(sptree,taxa2voucher_dic)
    path2_treeid_dic={}
    path_str_num_dic={}
    out=open('gd_loss_summary.txt','w')
    out.write('tree_ID'+'\t'+'gd_ID'+'\t'+'gd_support'+'\t'+'level'+'\t'+'species'+'\t'+'gene1'+'\t'+'gene2'+'\t'+'loss_path'+'\n')
    gd_id=1
    for tre_id,tre_path in tre_dic.items():
        t=PhyloTree(tre_path)
        t1=rename_input_tre(t, gene2new_named_gene_dic)
        path_str_num_lst,gd_id=get_path_str_with_count_num_lst(tre_id,gd_id,t1,renamed_sptree,out,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic,path2_treeid_dic)
        
        for i in path_str_num_lst :
            if i in path_str_num_dic:
                path_str_num_dic[i]+=1
            else:
                path_str_num_dic[i]=1
            
    out.close()
    return path_str_num_dic,path2_treeid_dic



if __name__ == "__main__":
    out='outfile'
    sptree=PhyloTree(sptree_path)
    num_sptree(sptree)
    tre_dic=read_and_return_dict(gf)

    os.makedirs(out, exist_ok=True)
    sp_dic,path2_treeid_dic=get_path_str_num_dic(tre_dic)

    split_dicts=split_dict_by_first_last_char(sp_dic)
    divide_path_results_into_individual_files_by_species(split_dicts,out)
    write_total_lost_path_counts_result(sp_dic)



