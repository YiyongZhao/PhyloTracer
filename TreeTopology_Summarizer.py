from __init__ import *
import os
  
def get_only_sps_tree(Phylo_t):
    for node in Phylo_t:
        node.name=node.name.split('_')[0]
    return Phylo_t
    
def get_species_counts_dic(node):
    species_repeats = {}
    for leaf in node:
        species = leaf.name
        if species in species_repeats:
            species_repeats[species] += 1
        else:
            species_repeats[species] = 1
    return species_repeats
    
def folding_tree(Phylo_t):  
    Phylo_t_c=Phylo_t.copy()
    for i in Phylo_t_c.traverse():
        i.add_feature('label',1)
    
    sp_set=set()
    for node in Phylo_t_c.traverse():
        if not node.is_leaf():
            sps=get_species_list(node)
            if len(set(sps))==1:
                node.name=sps[0]#+'_'+str(len(sps))
                #node.add_feature('label',str(len(sps)))
                node.label=len(sps)
                for child in node.get_children():
                    child.detach() 
            else:
                sp=tuple(get_species_set(node))
                if sp in sp_set:
                    continue
                else:
                    sp_set.add(sp)
                    if len(node.children) == 2:
                        left_child = node.children[0]
                        right_child = node.children[1]
                        node_up=node.up
                        if has_dup_clade(node):
                            s=get_species_counts_dic(node)
                            s1=str(list(set(s.values()))[0])
                            if len(left_child) >len(right_child):
                                right_child.label=int(s1)
                                if node_up:
                                    node.detach()
                                    node_up.add_child(right_child)
                                else:
                                    left_child.detach()
                            else:
                                left_child.label=int(s1)
                                if node_up:
                                    node.detach()
                                    node_up.add_child(left_child)
                                else:
                                    right_child.detach()      
    return Phylo_t_c

def get_species_set(node):
    return set(leaf.name for leaf in node.iter_leaves())
    
def has_dup_clade(tree):
    def compare_branch(node):
        left_child = node.children[0]
        right_child = node.children[1]
        if get_species_set(left_child) ==get_species_set(right_child):
        
            return True

    return compare_branch(tree)





def get_relative_clade_dic(relative_clade_dic,Phylo_t):
    for node in Phylo_t.traverse('postorder'):
        if not node.is_leaf():
            clade=node.write(format=9)
            num=get_num(node)
            if clade in relative_clade_dic :
                relative_clade_dic[clade]+=num
            else:
                relative_clade_dic[clade]=num
                
def statistical_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic):
    relative_clade_dic={}
    obtain_clade_dic={}
    for tre_path in tre_dic.values():
        Phylo_t0 = read_tree(tre_path)
        Phylo_t0 =root_tre_with_midpoint_outgroup(Phylo_t0)
        Phylo_t1=rename_input_tre(Phylo_t0,gene2new_named_gene_dic)
        Phylo_t1=get_only_sps_tree(Phylo_t0)
        Phylo_t2=folding_tree(Phylo_t1)
        get_relative_clade_dic(relative_clade_dic,Phylo_t2)
        get_obtain_clade_dic(obtain_clade_dic,Phylo_t2)
    #relative_new_dic=rejust_clade_dic(relative_clade_dic)
    #obtain_new_dic=rejust_clade_dic(obtain_clade_dic)
    relative_sorted_dict = dict(sorted(relative_clade_dic.items(), key=lambda x: len(x[0]), reverse=True))
    obtain_sorted_dict = dict(sorted(obtain_clade_dic.items(), key=lambda x: len(x[0]), reverse=True))
    with open ('Relative_Statistical_calde.txt','w') as f :
        for k,v in relative_sorted_dict.items():
            t=Tree(k)
            rename_input_tre(t,voucher2taxa_dic)
            f.write(t.write(format=9)+'\t'+str(v)+'\n')
    with open ('Obtain_Statistical_calde.txt','w') as f :
        for k,v in obtain_sorted_dict.items():
            t=Tree(k)
            rename_input_tre(t,voucher2taxa_dic)
            f.write(t.write(format=9)+'\t'+str(v)+'\n')


if __name__ == "__main__":
    tre_dic=read_and_return_dict('GF_list.txt')   
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    statistical_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic)
