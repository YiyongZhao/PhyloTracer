
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

def collapse_expand(Phylo_t:object,support_value:int)->object:
	for node in Phylo_t.traverse():
		if not node.is_leaf():
			if node.support <support_value:
			     node.delete()

	return Phylo_t

def collapse_expand_main(tre_dic:dict,support_value:int):
    dir_path = os.path.join(os.getcwd(), "output/collapse_expand_tree/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    for  k,v in tre_dic.items():
        t=Tree(v)
        t1=collapse_expand(t,support_value)
        tree_str=t1.write(format=0)
        write_tree_to_newick(tree_str,k,dir_path)



if __name__ == "__main__":
    
    tre_dic=read_and_return_dict('GF.txt')
    support_value=50
    collapse_expand_main(tre_dic,support_value,keep_sc,decimal_places)
    
