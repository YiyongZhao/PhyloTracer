
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

def collapse_nodes(Phylo_t: object, node_support_value: int) -> object:
    Phylo_t_c = Phylo_t.copy()
    
    for node in Phylo_t_c.traverse():
        if not node.is_leaf(): 
            if node.support < node_support_value:
                node.delete(preserve_branch_length=True)
    
    return Phylo_t_c


def collapse_expand_main(tre_dic:dict,support_value:int,revert:bool=False):
    dir_path = os.path.join(os.getcwd(), "collapse_expand_tree/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID,tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t=Tree(tre_path)
        t1=collapse_nodes(t,support_value)
        if revert:
            t1.resolve_polytomy(recursive=True)
            t1.sort_descendants("support")
        tree_str=t1.write(format=0)
        write_tree_to_newick(tree_str,tre_ID,dir_path)
        pbar.update(1)
    pbar.close()


if __name__ == "__main__":
    
    tre_dic=read_and_return_dict('GF.txt')
    support_value=50
    collapse_expand_main(tre_dic,support_value,keep_sc,decimal_places)
    
