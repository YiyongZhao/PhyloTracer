
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

def collapse_expand(Phylo_t:object,support_value:int)->object:
	for node in Phylo_t.tranverse():
		if not node.is_leaf():
			if node.support >support_value:
			     node.delete()

	return Phylo_t

if __name__ == "__main__":
    dir_path1 = os.path.join(os.getcwd(), "pruned_tree")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    tre_dic=read_and_return_dict('100_nosingle_GF_list.txt')
    support_value=50
    for  k,v in tre_dic.items():
    	t=Tree(v)
    	t1=collapse_expand(t)
    	write_tree_to_newick(t1,k,dir_path1)
