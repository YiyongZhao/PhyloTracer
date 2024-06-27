from __init__ import *

def trans_branch_length(Phylo_t: object, decimal_places: int=10) -> str:
    tree_str = Phylo_t.write(format=0, dist_formatter='%.{}f'.format(decimal_places))
    return tree_str

def write_tree_to_newick(tree_str:str,tre_ID:str,dir_path:str):
    with open(os.path.join(dir_path, tre_ID + '.nwk'),'w') as f :
        f.write(tree_str+'\n')

def branch_length_numeric_converter_main(tre_dic,decimal_places):
    dir_path = os.path.join(os.getcwd(), "converter_tree/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)
    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID,tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t=Tree(tre_path)
        if decimal_places==None:
            tree_str=trans_branch_length(t)
        else:
            tree_str=trans_branch_length(t,decimal_places)
        write_tree_to_newick(tree_str,tre_ID,dir_path)
        pbar.update(1)
    pbar.close()

if __name__ == "__main__":
    tre_dic=read_and_return_dict('test.txt') 
    branch_length_numeric_converter_main(tre_dic)
    



