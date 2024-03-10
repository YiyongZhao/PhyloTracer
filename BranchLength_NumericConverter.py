from __init__ import *

def trans_branch_length(Phylo_t:object)->str:
    tree_str=Phylo_t.write(format=0,dist_formatter='%.10f')
    return tree_str

def write_tree_to_newick(Phylo_t:object,tre_ID,dir_path1):
    tree_str=trans_branch_length(Phylo_t)
    with open(os.path.join(dir_path1, tre_ID + '.nwk'),'w') as f :
        f.write(tree_str+'\n')

def branch_length_numeric_converter_main(tre_dic):
    dir_path1 = os.path.join(os.getcwd(), "converter_tree")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    for tre_ID,tre_path in tre_dic.items():
        t=Tree(tre_path)
        write_tree_to_newick(t,tre_ID,dir_path1)
    

if __name__ == "__main__":
    tre_dic=read_and_return_dict('test.txt') 
    branch_length_numeric_converter_main(tre_dic)
    



