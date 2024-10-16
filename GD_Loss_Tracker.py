from __init__ import *
import re
def find_dup_node(Phylo_t:object)->list:#After searching for duplication events, the list of node names where duplication events occurred is as follows:
    events = Phylo_t.get_descendant_evol_events()
    dup_node = []
    for ev in events:
        if ev.etype == "D":
            i = ",".join(ev.in_seqs) + ',' + ",".join(ev.out_seqs)
            events_node_name_list = i.split(',')
            common_ancestor_node= Phylo_t.get_common_ancestor(events_node_name_list)
            dup_node.append(common_ancestor_node)
    return dup_node
    
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


def get_path_str_with_count_num_lst(tre_id,gd_id,genetree,sptree,out,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic,path2_treeid_dic):#统计一颗树中所有gd下不同的丢失情况
    renamed_sptree=rename_input_tre(sptree,taxa2voucher_dic)
    path_str_with_count_num_lst = []
    dup_node_name_list = find_dup_node(genetree)
    for i in dup_node_name_list:
        sp = get_species_set(i)
        clade_up = i.get_children()[0]
        clade_down = i.get_children()[1]
        up_tips='-'.join(clade_up.get_leaf_names())
        down_tips='-'.join(clade_down.get_leaf_names())

        if len(sp) != 1:
            max_clade2sp = renamed_sptree.get_common_ancestor(sp)
            sps_loss=set(max_clade2sp.get_leaf_names())-sp
            dup_sps=get_species_set(clade_up)&get_species_set(clade_down)
            gene_loss=get_gene_loss_list(clade_up,clade_down,dup_sps)

            dic,keys_with_zero_value=get_maptree_internal_node_name_count_dic(i,max_clade2sp,sptree)


            path_str_lst = get_tips_to_clade_path_lst(max_clade2sp,dic)
            for i in path_str_lst:
                if i in path2_treeid_dic:
                    path2_treeid_dic[i].append(f'{tre_id}-{gd_id}')
                else:
                    path2_treeid_dic[i]=[f'{tre_id}-{gd_id}']

            path_str_with_count_num_lst+=path_str_lst
            out.write(tre_id+'\t'+str(gd_id)+'\t'+max_clade2sp.name+'\t'+up_tips+'\t'+down_tips+'\t'+'-'.join(keys_with_zero_value)+'\t'+'-'.join(sorted(sps_loss))+'\t'+'-'.join(sorted(gene_loss))+'\n')
        else:
            out.write(tre_id+'\t'+str(gd_id)+'\t'+voucher2taxa_dic.get(get_maptree_name(renamed_sptree,sp),get_maptree_name(renamed_sptree,sp))+'\t'+up_tips+'\t'+down_tips+'\n')
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
    path2_treeid_dic={}
    path_str_num_dic={}
    out=open('gd_summary.txt','w')
    out.write(f'Tree_id\tGD_id\tSpecies_level\tGD_subclade1\tGD_subclade2\tNode_loss(twice)\tSpecies_Loss\tGene_Loss\n')
    gd_id=1
    for tre_id,tre_path in tre_dic.items():
        t=PhyloTree(tre_path)
        t1=rename_input_tre(t, gene2new_named_gene_dic)
        path_str_num_lst,gd_id=get_path_str_with_count_num_lst(tre_id,gd_id,t1,sptree,out,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic,path2_treeid_dic)
        for i in path_str_num_lst :
            if i in path_str_num_dic:
                path_str_num_dic[i]+=1
            else:
                path_str_num_dic[i]=1
            
    out.close()
    return path_str_num_dic,path2_treeid_dic

def write_total_lost_path_counts_result(sp_dic, path_dic):
    sorted_sp_keys = sorted(sp_dic.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_sp_dict = {k: sp_dic[k] for k in sorted_sp_keys}

    sorted_path_keys = sorted(path_dic.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_path_dict = {k: path_dic[k] for k in sorted_path_keys}

    with open('gd_loss_count_summary.txt', 'w') as f:
        f.write('GD Loss path\tGF count\n')
        processed_sp = set()
        for k, v in sorted_sp_dict.items():
            last_char = k.split('->')[-1].split('(')[0]
            if last_char not in processed_sp:
                f.write(f'\n{k}\t{v}\n')
                processed_sp.add(last_char)
            else:
                f.write(f'{k}\t{v}\n')

    with open('gd_loss_gf_count_summary.txt', 'w') as f1:
        f1.write('GD Loss path\tGF count\n')
        processed_path = set()
        for k1, v1 in sorted_path_dict.items():
            last_char1 = k1.split('->')[-1].split('(')[0]
            if last_char1 not in processed_path:
                f1.write(f'\n{k1}\t' + '\t'.join(map(str, v1)) + '\n')
                processed_path.add(last_char1)
            else:
                f1.write(f'{k1}\t' + '\t'.join(map(str, v1)) + '\n')

def proecee_start_node(file, sptree):
    with open(file, 'r') as f:
        species_list = [line.strip() for line in f.readlines()]
    try:
        common_ancestor_node = sptree.get_common_ancestor(species_list)
        return common_ancestor_node.name 
    except:
        print(f"Error: Unable to find a common ancestor for species: {species_list}")
        return None

def sort_dict_by_keys(input_dict):
    sorted_keys = sorted(input_dict.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_dict = {k: input_dict[k] for k in sorted_keys}
    return dict(sorted(sorted_dict.items(), key=lambda x: [int(n) for n in re.findall(r'\((\d+)\)', x[0])], reverse=True))

def write_gd_loss_to_file(filename, data_dict, key_filter, filter_value):
    with open(filename, 'w') as f:
        f.write('GD Loss path\tGF count\n')
        for k, v in data_dict.items():
            if key_filter(k) == filter_value:
                f.write(f'\n{k}\t{v}\n')

def write_gd_loss_info_of_strart_node(sp_dic, start_node, path_dic):
    sorted_dict1 = sort_dict_by_keys(sp_dic)
    sorted_dict2 = sort_dict_by_keys(path_dic)

    write_gd_loss_to_file('gd_loss_count_summary.txt', sorted_dict1, 
                           lambda k: k.split('->')[0].split('(')[0], start_node)
    
    write_gd_loss_to_file('gd_loss_gf_count_summary.txt', sorted_dict2, 
                           lambda k: k.split('->')[0].split('(')[0], start_node)

def write_gd_loss_info_of_species(sp_dic, species, path_dic):
    sorted_sp_dict = sort_dict_by_keys(sp_dic)
    sorted_path_dict = sort_dict_by_keys(path_dic)

    write_gd_loss_to_file('gd_loss_count_summary.txt', sorted_sp_dict, 
                           lambda k: k.split('->')[-1].split('(')[0], species)

    write_gd_loss_to_file('gd_loss_gf_count_summary.txt', sorted_path_dict, 
                           lambda k: k.split('->')[-1].split('(')[0], species)

def write_gd_loss_info_of_node_to_species(sp_dic, start_node, species, path_dic):
    sorted_dict1 = sort_dict_by_keys(sp_dic)
    sorted_dict2 = sort_dict_by_keys(path_dic)

    write_gd_loss_to_file('gd_loss_count_summary.txt', sorted_dict1, 
                           lambda k: (k.split('->')[0].split('(')[0], 
                                       k.split('->')[-1].split('(')[0]), 
                           (start_node, species))
    
    write_gd_loss_to_file('gd_loss_gf_count_summary.txt', sorted_dict2, 
                           lambda k: (k.split('->')[0].split('(')[0], 
                                       k.split('->')[-1].split('(')[0]), 
                           (start_node, species))


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



