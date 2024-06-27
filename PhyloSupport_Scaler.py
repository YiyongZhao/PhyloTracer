from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

def scale_support(Phylo_t: object, scale: str) -> object:
    max_support = max(node.support for node in Phylo_t.traverse())
    
    if max_support < 1:
        if scale == '100':
            for node in Phylo_t.traverse():
                node.support *= 100
        else:
            print('No need to scale down support, as all tree node supports are already < 1.')
    else:
        if scale != '100':
            for node in Phylo_t.traverse():
                node.support /= 100
        else:
            print('No need to scale support, as tree node supports are already in the range [1, 100].')
    
    return Phylo_t


def support_scaler_main(tre_dic,scale):
    dir_path = os.path.join(os.getcwd(), "support_scaler_tree/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID,tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t=Tree(tre_path)
        t1=scale_support(t,scale)
        tree_str=t1.write(format=0)
        write_tree_to_newick(tree_str,tre_ID,dir_path)
        pbar.update(1)
    pbar.close()

if __name__ == "__main__":
    scale=1
    tre_dic=read_and_return_dict('100_nosingle_GF_list.txt')
    support_scaler_main(tre_dic,scale)
