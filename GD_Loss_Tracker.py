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

def get_path_str_with_count_num_lst(tre_id,gd_id,genetree,sptree,out,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic):#统计一颗树中所有gd下不同的丢失情况
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

            out.write(tre_id+'\t'+str(gd_id)+'\t'+max_clade2sp.name+'\t'+up_tips+'\t'+down_tips+'\t'+'-'.join(sorted(sps_loss))+'\t'+'-'.join(sorted(gene_loss))+'\n')

            clade_up_set = get_maptree_internal_node_name_set(clade_up, max_clade2sp)
            clade_down_set = get_maptree_internal_node_name_set(clade_down, max_clade2sp)
            clade_up_lst = add_extra_maptree_node_name(clade_up_set, max_clade2sp)
            clade_down_lst = add_extra_maptree_node_name(clade_down_set, max_clade2sp)
            sp_list = list(clade_up_lst) + list(clade_down_lst)
            
            dic = get_maptree_node_count_dic(sp_list, max_clade2sp)

            path_str_lst = get_tips_to_clade_path_lst(max_clade2sp)
            for j in path_str_lst:
                s = '->'.join([f'{voucher2taxa_dic.get(key,key)}({dic[key]})' for key in j.split('->')])
                path_str_with_count_num_lst.append(s)
        else:
            out.write(tre_id+'\t'+str(gd_id)+'\t'+voucher2taxa_dic.get(get_maptree_name(renamed_sptree,sp),get_maptree_name(renamed_sptree,sp))+'\t'+up_tips+'\t'+down_tips+'\n')
        gd_id+=1

    set_path_str_with_count_num_lst = path_str_with_count_num_lst.copy()

    # temp_dic={}
    # for path in set_path_str_with_count_num_lst:
    #     if path.split('->')[-1].split('(')[0] in temp_dic:
    #         temp_dic[path.split('->')[-1].split('(')[0]].append(path)
    #     else:
    #         temp_dic[path.split('->')[-1].split('(')[0]]=[path]

    # newlist=[]
    # for k,v in temp_dic.items():
    #     templist=[]
    #     longest_element = max(v, key=len)
    #     if len(v)>1:
    #         for i in v:
    #             if i==longest_element:
    #                 templist.append(i)
    #                 continue
    #             copystr=longest_element[:len(longest_element)-len(i)]
    #             firstnum=int(i.split("->")[0].split("(")[-1].strip(")"))
    #             parts = copystr.split("->")
    #             lastnum = int(copystr.split("->")[-2].split("(")[-1].strip(")") if len(copystr.split("->")) > 1 else '0')
    #             if  firstnum>lastnum:
    #                 proecessstr=copystr[:-4]+str(firstnum)+')->'
    #                 newstr=proecessstr+i
    #             else:
    #                 newstr=copystr+i
    #             templist.append(newstr)
    #     else:
    #         templist.append(v[0])

        #newlist+=result
    result = sorted(set_path_str_with_count_num_lst, reverse=True)
        
    return result

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
    out.write(f'Tree_id\tGD_id\tSpecies_level\tGD_subclade1\tGD_subclade2\tspecies_Loss\tGene_Loss\n')
    gd_id=1
    for tre_id,tre_path in tre_dic.items():
        t=PhyloTree(tre_path)
        t1=rename_input_tre(t, gene2new_named_gene_dic)
        path_str_num_lst=get_path_str_with_count_num_lst(tre_id,gd_id,t1,sptree,out,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic,taxa2voucher_dic)
        for i in path_str_num_lst :
            if i in path_str_num_dic:
                path_str_num_dic[i]+=1
            else:
                path_str_num_dic[i]=1
            if i in path2_treeid_dic:
                path2_treeid_dic[i].append(tre_id)
            else:
                path2_treeid_dic[i]=[tre_id]
    out.close()
    return path_str_num_dic,path2_treeid_dic

def divide_path_results_into_individual_files_by_species(split_dicts, out_dir):
    for file_name,dic in split_dicts.items():  
        with open(f"{out_dir}/{file_name}.txt", 'w') as output_file:
            for path, num in dic.items():
                output_file.write(f"{path}\t{num}\n")

def write_total_lost_path_counts_result(sp_dic):
    sorted_keys = sorted(sp_dic.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_dict = {k: sp_dic[k] for k in sorted_keys}
    with open('gd_loss_count_summary.txt','w') as f :
        f.write(f'GD Loss path\tGF count')
        precess=set()
        for k,v in sorted_dict.items():
            last_char = k.split('->')[-1].split('(')[0]
            if last_char not in precess:
                f.write('\n')
                f.write(k+'\t'+str(v)+'\n')
                precess.add(last_char)
            else:
                f.write(k+'\t'+str(v)+'\n')

def write_total_lost_path_treeid_result(sp_dic):
    sorted_keys = sorted(sp_dic.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_dict = {k: sp_dic[k] for k in sorted_keys}
    with open('gd_loss_gf_count_summary.txt','w') as f :
        f.write(f'GD Loss path\tGF count')
        precess=set()
        for k,v in sorted_dict.items():
            last_char = k.split('->')[-1].split('(')[0]
            if last_char not in precess:
                f.write('\n')
                f.write(k+'\t'+'\t'.join(v)+'\n')
                precess.add(last_char)
            else:
                f.write(k+'\t'+'\t'.join(v)+'\n')

def proecee_start_node(file, sptree):
    with open(file, 'r') as f:
        species_list = [line.strip() for line in f.readlines()]
    try:
        common_ancestor_node = sptree.get_common_ancestor(species_list)
        return common_ancestor_node.name 
    except:
        print(f"Error: Unable to find a common ancestor for species: {species_list}")
        return None

def write_gd_loss_info_of_strart_node(sp_dic,start_node):
    sorted_keys = sorted(sp_dic.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_dict = {k: sp_dic[k] for k in sorted_keys}
    with open('gd_loss_count_summary.txt','w') as f :
        f.write(f'GD Loss path\tGF count')
        for k,v in sorted_dict.items():
            first_char = k.split('->')[0].split('(')[0]
            if first_char==start_node:
                f.write('\n')
                if isinstance(v, (list, tuple)):
                    f.write(k + '\t' + '\t'.join(v) + '\n')
                else:
                    f.write(k + '\t' + str(v) + '\n')  
                    
def write_gd_loss_info_of_species(sp_dic,species):
    sorted_keys = sorted(sp_dic.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_dict = {k: sp_dic[k] for k in sorted_keys}
    with open('gd_loss_count_summary.txt','w') as f :
        f.write(f'GD Loss path\tGF count')
        for k,v in sorted_dict.items():
            first_char = k.split('->')[-1].split('(')[0]
            if first_char==species:
                f.write('\n')
                if isinstance(v, (list, tuple)):
                    f.write(k + '\t' + '\t'.join(v) + '\n')
                else:
                    f.write(k + '\t' + str(v) + '\n')  



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



