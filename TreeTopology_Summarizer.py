from __init__ import *

def get_only_sps_tree(Phylo_t):
    Phylo_t_c=Phylo_t.copy()
    for node in Phylo_t_c:
        node.name=node.name.split('_')[0]
    return Phylo_t_c
    
def get_max_tree(trees):
    max_leaf_count = 0
    max_leaf_count_tree = None

    for i in trees:
        leaf_count = len(i.get_leaves())
        if leaf_count > max_leaf_count:
            max_leaf_count = leaf_count
            max_leaf_count_tree = i

    return max_leaf_count_tree

def process(trees, result_dict):
    if len(trees) == 1:
        result_dict[trees[0].write(format=9)] = 1
        return []
    else:
        same_rf_trees = []
        different_trees = []
        first_tree = trees[0]
        for i in trees:
                
            rf = first_tree.robinson_foulds(i)[0]
            if rf == 0:
                same_rf_trees.append(i)
            else:
                different_trees.append(i)
           
            
        #max_tree = get_max_tree(same_rf_trees)
        #result_dict[max_tree.write(format=9)] = len(same_rf_trees)
        result_dict[max_tree.write(format=9)] = same_rf_trees
            
        #o.write(max_tree.write(format=9)+'\t'+str(len(same_rf_trees))+'\n')
        return different_trees
      
def process_tree(trees,result_dict):
    remaining_trees = trees
    while len(remaining_trees) >= 1:
        remaining_trees = process(remaining_trees, result_dict)
    
    
if __name__ == "__main__":
    tre_dic=read_and_return_dict('GF_list.txt')   
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    statistical_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic)
