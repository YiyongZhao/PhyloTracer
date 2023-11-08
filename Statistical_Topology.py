from __init__ import *
import os

def get_empty_count_dict(sptree: object) -> dict:
    empty_count_dic={}
    for node in sptree.traverse():
        if not node.is_leaf():
            clade=node.write(format=9)
            empty_count_dic[clade]=0
    return empty_count_dic

def switch_tree(phylo_tree):
    
    down_child = phylo_tree.get_children()[1]
    up_child = phylo_tree.get_children()[0]
    
    new_tree = Tree()
    new_tree.add_child(down_child)
    new_tree.add_child(up_child)
    
    return new_tree

def get_only_sps_tree(Phylo_t):
    for node in Phylo_t:
        node.name=node.name.split('_')[0]
    return Phylo_t

def folding_tree(Phylo_t):
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            sps=get_species_list(node)
            if len(set(sps))==1:
                node.name=sps[0]
                for child in node.get_children():
                    child.detach() 
    return Phylo_t

def statistical_clade(clade_dic,Phylo_t):
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            clade=node.write(format=9)
            if clade in clade_dic :
                clade_dic[clade]+=1
            else:
                clade_dic[clade]=1

def get_multiplier(t):
    multiplier=1
    for i in t.traverse():
        if not i.is_leaf():
            sps=get_species_list(i)
            if len(set(sps))==1:
                 multiplier=multiplier*len(sps)
    return multiplier
                
def rejust_clade_dic(clade_dic):
    new_dic={}
    for k,v in clade_dic.items():
        t=Tree(k)
        t1=switch_tree(t)
        t2=t1.write(format=9)
        if k and t2 not in new_dic :
            new_dic[k]=v
        if t2 in new_dic :
            new_dic[t2]=new_dic[t2]+v
            
    return new_dic


def statistical_main():
    relative_clade_dic={}
    obtain_clade_dic={}
    for tre_path in tre_dic.values():
        Phylo_t0 = read_tree(tre_path)
        Phylo_t0 =root_tre_with_midpoint_outgroup(Phylo_t0)
        Phylo_t1=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
        Phylo_t1=get_only_sps_tree(Phylo_t0)
        Phylo_t2=folding_tree(Phylo_t1)
        statistical_clade(relative_clade_dic,Phylo_t2)
        statistical_clade(obtain_clade_dic,Phylo_t1)
    relative_new_dic=rejust_clade_dic(relative_clade_dic)
    obtain_new_dic=rejust_clade_dic(obtain_clade_dic)
    with open ('Relative_Statistical_calde.txt','w') as f :
        for k,v in relative_new_dic.items():
            t=Tree(k)
            rename_input_tre(t,voucher2taxa_dic)
            f.write(t.write(format=9)+'\t'+str(v)+'\n')
    with open ('Obtain_Statistical_calde.txt','w') as f :
        for k,v in obtain_new_dic.items():
            t=Tree(k)
            if len(set(get_species_list(t))) !=1:
                s=get_multiplier(t)
                folding_tree(t)
                rename_input_tre(t,voucher2taxa_dic)
                f.write(t.write(format=9)+'\t'+str(v*s)+'\n')


if __name__ == "__main__":
  tre_dic=read_and_return_dict('GF_list.txt')   
  gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
  statistical_main()
