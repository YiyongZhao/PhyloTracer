from __init__ import *
from Tree_Visualizer import *

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

    


def write_gene_duplication_events(sp_dic,filename, tre_dic, support,dup_species_percent, dup_species_num, sptree,gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic):
    with open(filename,'w') as file:
        file.write('tree_ID'+'\t'+'gd_id'+'\t'+'gd_support'+'\t'+'gene1'+'\t'+'gene2'+'\t'+'level'+'\t'+'species'+'\t'+'GD_dup_sps'+'\t'+'dup_ratio'+'\t'+'gd_type'+'\t'+'comment'+'\n') 
        gd_num=1
        
        for tre_ID, tre_path in tre_dic.items():
            Phylo_t0 = read_phylo_tree(tre_path)
            Phylo_t0=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
            #Phylo_t1=root_tre_with_midpoint_outgroup(Phylo_t0)
            num_tre_node(Phylo_t0)
            dup_node_name_list = find_dup_node(Phylo_t0)
            #Phylo_t0.write(outfile='num_tree/' + str(tre_ID) + '.nwk', format=1)
            tree_dic={}
            for i in dup_node_name_list:
                
                clade=Phylo_t0&i
                sps=get_species_list(clade)
                parent=sptree.get_common_ancestor(sps)
                dup_sps=get_duplicated_species_num(clade)
                dup_percent=dup_sps/len(set(sps))
                model=get_model(clade,sptree)
                gene_pair1=gene_pair(clade)
                null='null'
                tree_dic[parent.name]=model
                
                for j in gene_pair1:
                    file.write(tre_ID+'\t'+str(gd_num)+'\t')
                    a,b=j.split('-')
                    if a.split('_')[0] ==b.split('_')[0]:
                        c=a.split('_')[0]
                    else:
                        if a=='null':
                            c=b.split('_')[0]
                        else:
                            c=a.split('_')[0]
        
                    file.write(str(clade.support)+'\t'+new_named_gene2gene_dic.get(a, null)+'\t'+new_named_gene2gene_dic.get(b, null)+'\t'+parent.name+'\t'+voucher2taxa_dic[c]+'\t'+str(dup_sps)+'\t'+str(round(dup_percent,2))+'\t'+model+'\t'+'-'+'\t''\n')
                gd_num+=1
            sp_dic.append(tree_dic)
 
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
    
def batch_gfs_traverse(tre_dic: dict, support_value: int, empty_count_dic: dict,sptree:object,gene2new_named_gene_dic) -> dict:
    for tre_path in tre_dic.values():
        Phylo_t0 = read_phylo_tree(tre_path)
        Phylo_t1=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
        #Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
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

def mark_sptree(sptree:object,empty_count_dic:dict,voucher2taxa_dic)->object:
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
    sptree=rename_input_tre(sptree,voucher2taxa_dic)
    return sptree.render('species_tree_GD_Detector.PDF')
    
def rename_species_tree(sptree:object, voucher2taxa_dic:dict)->object:
    for node in sptree:
        key_to_find = node.name
        if key_to_find in voucher2taxa_dic.values():
            new_name = next(key for key, value in voucher2taxa_dic.items() if value == key_to_find)
            node.name = new_name           
    return sptree
    
def num_tre_node(Phylo_t:object)->object:#Numbering the nodes in the tree.
    i = 1
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            node.name = "N" + str(i)
            i += 1
    return Phylo_t

def get_species_list(Phylo_t:object)->list:
    leaves = Phylo_t.get_leaf_names()
    species_lst = []
    for leaf in leaves:
        species = leaf.split("_")[0]  # Assuming that the species name is located in the first part of the node name (leaf) with “_” as the delimiter
        species_lst.append(species)
    return species_lst


def gene_pair(clade):
    result_pairs = set()

    child1, child2 = clade.get_children()[0].get_leaf_names(), clade.get_children()[1].get_leaf_names()

    for tip1 in child1:
        matching_tips = [tip2 for tip2 in child2 if tip1.split('_')[0] == tip2.split('_')[0]]
        
        if matching_tips:
            result_pairs.update(f"{tip1}-{tip2}" for tip2 in matching_tips)
        else:
            result_pairs.add(f"{tip1}-null")

    
    for tip2 in child2:
        if all(tip2.split('_')[0] != tip1.split('_')[0] for tip1 in child1):
            result_pairs.add(f"null-{tip2}")

    return result_pairs

def find_tre_dup(Phylo_t:object) -> list:   #seperator either "@" or "_"
    tre_ParaL=[]
    GF_leaves_S = set(Phylo_t.get_leaf_names())
    events = Phylo_t.get_descendant_evol_events()
    for ev in events:
        if ev.etype == "D":
            tre_ParaL.append(",".join(ev.in_seqs)+"<=>"+",".join(ev.out_seqs))
    return tre_ParaL,GF_leaves_S    


def count_values_in_list_of_dicts(lst_of_dicts):
    result_dict = {}

    for i in lst_of_dicts:
        for k, v in i.items():
            if k not in result_dict:
                result_dict[k] = {v: 1}
            else:
                if v in result_dict[k]:
                    result_dict[k][v] += 1
                else:
                    result_dict[k][v] = 1

    return result_dict

def save_matrix_to_csv(result_dict, filename='output_matrix.csv'):

    df = pd.DataFrame(result_dict).fillna(0).astype(int)
    df = df.transpose()

    df.to_csv(filename, index_label='')
    

if __name__ == "__main__":
    support=50
    dup_species_percent = 0.5
    dup_species_num = 2
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap.txt")
    sptree=PhyloTree('30sptree.nwk')
    sptree=rename_species_tree(sptree, voucher2taxa_dic)
    num_tre_node(sptree)
    tre_dic=read_and_return_dict('GF.txt')
    filename = 'result.txt'
    sp_dic=[]
    write_gene_duplication_events(sp_dic,filename, tre_dic, support,dup_species_percent, dup_species_num,sptree,gene2new_named_gene_dic,new_named_gene2gene_dic,voucher2taxa_dic)
    empty_count_dic=get_empty_count_dict(sptree)
    empty_count_dic=batch_gfs_traverse(tre_dic, support, empty_count_dic,sptree,gene2new_named_gene_dic) 
    mark_sptree(sptree,empty_count_dic,voucher2taxa_dic)
    p=count_values_in_list_of_dicts(sp_dic)    
    save_matrix_to_csv(p)
   
