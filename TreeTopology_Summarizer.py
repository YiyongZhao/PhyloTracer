from __init__ import *

def get_max_tree(trees):
    max_leaf_count = 0
    max_leaf_count_tree = None

    for i in trees:
        leaf_count = len(i.get_leaves())
        if leaf_count > max_leaf_count:
            max_leaf_count = leaf_count
            max_leaf_count_tree = i

    return max_leaf_count_tree

def process(trees, result_dict,dic):
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
           
            
        max_tree = get_max_tree(same_rf_trees)
        result_dict[max_tree.write(format=9)] = len(same_rf_trees)
        dic[max_tree.write(format=9)] = same_rf_trees
            
        #o.write(max_tree.write(format=9)+'\t'+str(len(same_rf_trees))+'\n')
        return different_trees
      
def process_tree(trees, output_file,dic):
    result_dict = {}
    remaining_trees = trees
    while len(remaining_trees) >= 1:
        remaining_trees = process(remaining_trees, result_dict,dic)
    sorted_result = sorted(result_dict.items(), key=lambda x: len(x[0]), reverse=True)
    with open(output_file, 'w') as file:
        for k, v in sorted_result:
            file.write(f"{k}\t{v}\n")
    
if __name__ == "__main__":
    tre_dic=read_and_return_dict('GF_list.txt')   
    gene2new_named_gene_dic, new_named_gene2gene_dic,voucher2taxa_dic=gene_id_transfer("imap")
    statistical_main(tre_dic,gene2new_named_gene_dic,voucher2taxa_dic)
