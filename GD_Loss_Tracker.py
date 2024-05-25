from __init__ import *

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
    
def get_maptree_internal_node_name_set(clade,sptree):
    names=set()
    for i in clade.traverse():
        if not i.is_leaf():
            sp_set=get_species_set(i)
            if len(sp_set) !=1 :
                sp=sptree.get_common_ancestor(sp_set)
                names.add(sp.name)
            else:
                sp_set=list(sp_set)[0]
                sp=sptree&sp_set
                names.add(sp.name)
        else:
            sp_set=i.name.split('_')[0]
            sp=sptree&sp_set
            names.add(sp.name)
    return names

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
    return dic 

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
    path_str='->'.join(re_nodes)
    return path_str

def get_tips_to_clade_path_lst(taget_node:object)->list:#taget_node也就是gd node
    path_str_lst=[]
    for i in taget_node:
        path_str=get_two_nodes_path_str(i,taget_node)
        path_str_lst.append(path_str)
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

def get_path_str_with_count_num_lst(tre_id,gd_id,genetree,sptree,out):#统计一颗树中所有gd下不同的丢失情况
    path_str_with_count_num_lst = []
    dup_node_name_list = find_dup_node(genetree)
    for i in dup_node_name_list:
        sp = get_species_set(i)
        clade_up = i.get_children()[0]
        clade_down = i.get_children()[1]
        up_tips='-'.join(clade_up.get_leaf_names())
        down_tips='-'.join(clade_down.get_leaf_names())

        if len(sp) != 1:
            max_clade2sp = sptree.get_common_ancestor(sp)
            out.write(tre_id+'\t'+str(gd_id)+'\t'+max_clade2sp.name+'\t'+up_tips+'\t'+down_tips+'\n')

            clade_up_set = get_maptree_internal_node_name_set(clade_up, max_clade2sp)
            clade_down_set = get_maptree_internal_node_name_set(clade_down, max_clade2sp)
            clade_up_lst = add_extra_maptree_node_name(clade_up_set, max_clade2sp)
            clade_down_lst = add_extra_maptree_node_name(clade_down_set, max_clade2sp)
            sp_list = list(clade_up_lst) + list(clade_down_lst)
            
            dic = get_maptree_node_count_dic(sp_list, max_clade2sp)

            path_str_lst = get_tips_to_clade_path_lst(max_clade2sp)
            for j in path_str_lst:
                s = '->'.join([f'{key}({dic[key]})' for key in j.split('->')])
                path_str_with_count_num_lst.append(s)
        else:
            out.write(tre_id+'\t'+str(gd_id)+'\t'+get_maptree_name(sptree,sp)+'\t'+up_tips+'\t'+down_tips+'\n')
        gd_id+=1
    return path_str_with_count_num_lst

def num_sptree(sptree):
    n=0
    for i in sptree.traverse('postorder'):
        if not i.is_leaf():
            i.name='S'+str(n)
            n+=1
    sptree.sort_descendants()
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

def get_path_str_num_dic(tre_dic,sptree):
    path_str_num_dic={}
    out=open('gd_summary.txt','w')
    gd_id=1
    for tre_id,tre_path in tre_dic.items():
        t=PhyloTree(tre_path)
        path_str_num_lst=get_path_str_with_count_num_lst(tre_id,gd_id,t,sptree,out)
        for i in path_str_num_lst :
            if i in path_str_num_dic:
                path_str_num_dic[i]+=1
            else:
                path_str_num_dic[i]=1
    out.close()
    return path_str_num_dic

def divide_path_results_into_individual_files_by_species(split_dicts, out_dir):
    for file_name,dic in split_dicts.items():  
        with open(f"{out_dir}/{file_name}.txt", 'w') as output_file:
            for path, num in dic.items():
                output_file.write(f"{path}\t{num}\n")

def write_total_lost_path_counts_result(sp_dic):
    with open('gd_loss_summary.txt','w') as f :
        for k,v in sp_dic.items():
            f.write(k+'\t'+str(v)+'\n')

      
if __name__ == "__main__":
    out='outfile'
    sptree=PhyloTree(sptree_path)
    num_sptree(sptree)
    tre_dic=read_and_return_dict(gf)

    os.makedirs(out, exist_ok=True)
    sp_dic=get_path_str_num_dic(tre_dic)
    split_dicts=split_dict_by_first_last_char(sp_dic)
    divide_path_results_into_individual_files_by_species(split_dicts,out)
    write_total_lost_path_counts_result(sp_dic)

    

