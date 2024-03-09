from __init__ import *

def get_scientific_to_float_dic(Phylo_t:object)->dict:
    scientific_to_float_dic={}
    for node in Phylo_t.traverse():
        if 'e' in str(node.dist):
            float_num="{:.10f}".format(node.dist)
            scientific_to_float_dic[str(node.dist)]=str(float_num)
    return scientific_to_float_dic

def trans_branch_length(Phylo_t:object)->str:
    tree_str=Phylo_t.write()
    scientific_to_float_dic=get_scientific_to_float_dic(Phylo_t)
    for scientific_str in scientific_to_float_dic :
        if scientific_str in tree_str :
            tree_str=tree_str.replace(scientific_str,scientific_to_float_dic[scientific_str])
    return tree_str

def write_tree_to_newick(Phylo_t:object,tre_ID):
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
        write_tree_to_newick(t,tre_ID)
    

if __name__ == "__main__":
    tre_dic=read_and_return_dict('GF_list.txt') 
    branch_length_numeric_converter_main(tre_dic)
    


