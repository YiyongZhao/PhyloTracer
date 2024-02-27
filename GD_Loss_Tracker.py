from __init__ import *

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

def prcess(lst,date):
    for i in lst:
        k, v = i.split('-')
        date[k][v]+=1
        
def get_mapp_list(clade,sptree):
    lst=set()
    for i in clade.traverse():
        if not i.is_leaf():
            sp_lst=get_species_set(i)
            if len(sp_lst) !=1 :
                sp=sptree.get_common_ancestor(sp_lst)
                lst.add(sp.name)
            else:
                sp_lst=list(sp_lst)[0]
                sp=sptree&sp_lst
                lst.add(sp.name)
                continue
        else:
            sp_lst=i.name.split('_')[0]
            sp=sptree&sp_lst
            lst.add(sp.name)
    return lst

def get_up_down_lst(max_clade,sptree):
    clade_up=max_clade.get_children()[0]
    clade_up_set=get_mapp_list(clade_up,sptree)
    clade_down=max_clade.get_children()[1]
    clade_down_set=get_mapp_list(clade_down,sptree)
    up_down_lst=list(clade_up_set)+list(clade_down_set)
    return up_down_lst 

def get_summary_dic(max_clade,max_clade2sp,sptree):
    up_down_lst=get_up_down_lst(max_clade,sptree)
    dic={i.name:0 for i in max_clade2sp.traverse()}
    for i in up_down_lst :
        if i in dic :
            dic[i]+=1
    return dic 

def sort_dic(dic):
    def custom_sort(key):
        if key.startswith('N'):
            return (0, key)
        else:
            return (1, key)
    sorted_dic={k: dic[k] for k in sorted(dic, key=custom_sort)}
    return sorted_dic
  
def main():
    for tre_ID, tre_path in tre_dic.items(): 
        Phylo_t0 = read_phylo_tree(tre_path)
        num_tre_node(Phylo_t0)
        dup_node_name_list = find_dup_node(Phylo_t0)
        if len(dup_node_name_list) ==0:
            continue
        i=dup_node_name_list[0]
        max_clade=Phylo_t0&i
        sp=get_species_set(max_clade)
        max_clade2sp=sptree.get_common_ancestor(sp)



        e=get_summary_dic(max_clade,max_clade2sp,sptree)
        f1=sort_dic(e)
        f2=[k+'-'+str(v) for k,v in f1.items()]
        #f3=[f2[0]]
        #prcess(f3,date)
        prcess(f2,date)
      
if __name__ == "__main__":
    tre_dic=read_and_return_dict('GF_list.txt')   
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    pattern_dic = {'2':0,'1':0,'0':0}
    sp_dic={i.name:pattern_dic for i in sptree.traverse()}
    model=['2','1','0']
    node=sort_dic(sp_dic)
    nodes=node.keys()
    date=pd.DataFrame(columns=nodes,index=model).fillna(0)
    main()
    date.to_csv('output.csv')
    

