from __init__ import *

def find_dup_node(Phylo_t:object)->list:#After searching for duplication events, the list of node names where duplication events occurred is as follows:
    events = Phylo_t.get_descendant_evol_events()
    dup_node_name_list = []
    for ev in events:
        if ev.etype == "D":
            i = ",".join(ev.in_seqs) + ',' + ",".join(ev.out_seqs)
            events_node_name_list = i.split(',')
            common_ancestor_node_name = Phylo_t.get_common_ancestor(events_node_name_list)
            dup_node_name_list.append(common_ancestor_node_name.name)
    return dup_node_name_list

def get_duplicated_species_num(node):
    species_lst=get_species_list(node)
    duplicates = []
    unique_elements = set()

    for element in species_lst:
        if element in unique_elements:
            duplicates.append(element)
        else:
            unique_elements.add(element)

    duplicate_count = len(set(duplicates))

    return duplicate_count

def get_intraspecies_and_interspecies_node_list(dup_node_name_list,Phylo_t0):
    intraspecies_node_list=[]
    interspecies_node_list=[]
    for node in dup_node_name_list:
        clade=Phylo_t0&node
        sps=get_species_list(clade)
        if len(set(sps))==1 :
            intraspecies_node_list.append(node)
        else:
            interspecies_node_list.append(node)
    return intraspecies_node_list,interspecies_node_list
    
def get_target_node_list(dup_node_name_list,Phylo_t0,dup_percent,dup_species_num):
    target_node_list=[]
    for node in dup_node_name_list:
        clade=Phylo_t0&node
        sps=get_species_list(clade)
        if get_duplicated_species_num(clade)/len(set(sps)) >=dup_percent and get_duplicated_species_num(clade)>=dup_species_num:
            target_node_list.append(node)
    return target_node_list

def write_gene_duplication_events(filename, tre_dic, support,dup_species_percent, dup_species_num, sptree):
    with open(filename,'w') as file:
        file.write('tre_ID'+'\t'+'gd_node_id'+'\t'+'level'+'\t'+'gd_support'+'\t'+'dup_sps'+'\t'+'dup_percent'+'\t'+'gd_type'+'\t'+'comment'+'\n') 
        for tre_ID, tre_path in tre_dic.items():
            
            Phylo_t0 = read_phylo_tree(tre_path)
            num_tre_node(Phylo_t0)
            dup_node_name_list = find_dup_node(Phylo_t0)
            #Phylo_t0.write(outfile='num_tree/' + str(tre_ID) + '.nwk', format=1)
            for i in dup_node_name_list:
                file.write(tre_ID+'\t'+i+'\t')
                clade=Phylo_t0&i
                sps=get_species_list(clade)
                dup_sps=get_duplicated_species_num(clade)
                dup_percent=dup_sps/len(set(sps))
                model=get_model(clade,sptree)
                if len(set(sps)) ==1:
                    file.write(sps[0]+'\t'+str(clade.support)+'\t'+str(dup_sps)+'\t'+str(round(dup_percent,2))+'\t'+model+'\t'+'Intraspecific_gene_duplication_event'+'\n')
                else:
                    if clade.support >=support and dup_sps>=dup_species_num and dup_percent>=dup_species_percent:
                        file.write(clade.name+'\t'+str(clade.support)+'\t'+str(dup_sps)+'\t'+str(round(dup_percent,2))+'\t'+model+'\t'+'Targetspecific_gene_duplication_event'+'\n')
                    else:
                        file.write(clade.name+'\t'+str(clade.support)+'\t'+str(dup_sps)+'\t'+str(round(dup_percent,2))+'\t'+model+'\t'+'Interspecific_gene_duplication_event'+'\n')
 
def get_model(clade,sptree):
    sps=get_species_list(clade)
    sps_clade=sptree.get_common_ancestor(set(sps))
    for leaf in sps_clade.get_children()[0]:
        leaf.add_feature('label','A')
    for leaf in sps_clade.get_children()[1]:
        leaf.add_feature('label','B')
    for j in clade:
        species=j.name.split('_')[0]
        clade1=sps_clade&species
        if clade1:
            j.add_feature('label',clade1.label)
    up_clade=''
    for j in clade.get_children()[0]:
        up_clade+=j.label
    up_clade=up_clade+'<=>' 
    for j in clade.get_children()[1]:
        up_clade+=j.label
    clade_up=set(up_clade.split('<=>')[0])
    clade_down=set(up_clade.split('<=>')[1])
    clade_model=''.join(clade_up)+'<=>'+''.join(clade_down)
    return clade_model

def get_model_dic(interspecies_node_list,genetree,sptree):
    model_dic={}
    for i in interspecies_node_list:
        clade=genetree&i
        s=get_model(clade,sptree)
        model_dic.setdefault(s, []).append(i)
    return model_dic

def get_empty_count_dict(sptree: object) -> dict:
    empty_count_dic = {node.name: 0 for node in sptree.traverse()}
    return empty_count_dic

def judge_support(support,support_value):
    if  support <=1 and 0.5 <=support_value <=1:
        if support>support_value:
            return True
        else:
            return False
        
    elif support <=1 and 50 <= support_value <=100:
        support_value=support_value/100
        if support>support_value:
            return True
        else:
            return False
    elif support > 1 and 0.5 <=support_value <=1:
        support_value=support_value*100
        if support>support_value:
            return True
        else:
            return False
    elif support > 1 and 50 <=support_value <=100:
        if support>support_value:
            return True
        else:
            return False
    
def batch_gfs_traverse(tre_dic: dict, support_value: int, empty_count_dic: dict,sptree:object) -> dict:
    for tre_path in tre_dic.values():
        Phylo_t0 = read_tree(tre_path)
        Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
        
        tre_ParaL,GF_leaves_S =find_tre_dup(Phylo_t1)
        for dup in tre_ParaL:
            leafs=dup.replace('<=>',',')
            leaf_list=leafs.split(',')
            ancestor=Phylo_t1.get_common_ancestor(leaf_list)
            if judge_support(ancestor.support,support_value):
                sps_list=set([leaf[0:3] for leaf in leaf_list])
                if len(sps_list)>1:
                    parent=sptree.get_common_ancestor(sps_list)
                    empty_count_dic[parent.name]+=1
                else:
                    single_sps=list(sps_list)[0]
                    parent=sptree&single_sps
                    empty_count_dic[parent.name]+=1
    return empty_count_dic

def mark_sptree(sptree:object,empty_count_dic:dict)->object:
    sptree.ladderize()
    sptree.sort_descendants("support")
    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        node.set_style(nstyle)
        if not node.is_leaf() and node.name in empty_count_dic:
            leafs=node.get_leaf_names()
            num_list=[empty_count_dic[leaf] for leaf in leafs]
            num_list.append(empty_count_dic[node.name])
            num=str(sum(num_list))
            node.add_face(TextFace(num+' GD', fsize=5, fgcolor="red"), column=0, position="branch-top")
            
    return sptree.render('species_tree_GD_Detector.PDF')
    
if __name__ == "__main__":
    sptree=PhyloTree('F:/a/96tree/7sp.nwk')
    num_tre_node(sptree)
    tre_dic=read_and_return_dict('F:/a/96tree/96.txt')
    empty_count_dic=get_empty_count_dict(sptree)
    empty_count_dic=batch_gfs_traverse(tre_dic, support_value, empty_count_dic,sptree) 
    mark_sptree(sptree,empty_count_dic)
    filename = '30sp_31131tree_GD.txt'
    support=50
    dup_species_percent = 0.5
    dup_species_num = 2
    write_gene_duplication_events(filename, tre_dic, support,dup_species_percent, dup_species_num,sptree)
