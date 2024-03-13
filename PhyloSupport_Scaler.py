
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick
def scale_support(Phylo_t:object,scale=True):
    max_support=max([i.support for i in Phylo_t.traverse()])
    if max_support<1:
        if scale:
            for node in Phylo_t.traverse():
                node.support=node.support*100
            return Phylo_t
                
        else:
            return Phylo_t
            print('there is no need to scale_down_support ,because the tree node support is all <1')
            
            
    else:
        if scale:
            
            
            return Phylo_t
            print('there is no need to scale_support ,because the tree node support is in range of [1,100]')
        else:
            for node in Phylo_t.traverse():
                node.support=node.support/100
            
            return Phylo_t

def support_scaler_main(tre_dic,scale):
    dir_path1 = os.path.join(os.getcwd(), "support_scaler_tree")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    for tre_ID,tre_path in tre_dic.items():
        t=Tree(tre_path)
        t1=scale_support(t,scale)
        write_tree_to_newick(t1,tre_ID,dir_path1)

if __name__ == "__main__":
    scale=True
    tre_dic=read_and_return_dict('100_nosingle_GF_list.txt')
    support_scaler_main(tre_dic,scale)
