
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

def scale_support(Phylo_t: object, scale: str = '100') -> object:
    """
    Scale the support values of all nodes in a phylogenetic tree.

    Args:
        Phylo_t (object): The phylogenetic tree object, must have 'traverse' and 'support' attributes.
        scale (str or int): Target scale, '100' to scale up to [1,100], otherwise scale down to [0,1].

    Returns:
        object: The tree object with scaled support values.

    Raises:
        ValueError: If scale is not a valid value.
    """
    if not hasattr(Phylo_t, "traverse"):
        raise AttributeError("Phylo_t must have a 'traverse' method.")
    try:
        max_support = max(getattr(node, 'support', 0) for node in Phylo_t.traverse())
    except Exception as e:
        raise RuntimeError(f"Failed to get support values: {e}")

    # 支持 scale 为字符串或数字
    if isinstance(scale, str) and scale.isdigit():
        scale = int(scale)

    if max_support < 1:
        if scale == 100:
            for node in Phylo_t.traverse():
                if hasattr(node, 'support'):
                    node.support *= 100
        else:
            print('No need to scale down support, as all tree node supports are already < 1.')
    else:
        if scale != 100:
            for node in Phylo_t.traverse():
                if hasattr(node, 'support'):
                    node.support /= 100
        else:
            print('No need to scale support, as tree node supports are already in the range [1, 100].')
    return Phylo_t

def support_scaler_main(tre_dic: dict, scale: str = '100') -> None:
    """
    Batch process trees to scale their support values and write to files.

    Args:
        tre_dic (dict): Dictionary mapping tree IDs to tree file paths.
        scale (str or int): Target scale for support values.

    Raises:
        ValueError: If tre_dic is not a dict.
        OSError: If directory or file operations fail.
    """
    if not isinstance(tre_dic, dict):
        raise ValueError("tre_dic must be a dictionary.")

    dir_path = os.path.join(os.getcwd(), "support_scaler_tree/")
    try:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path)
    except Exception as e:
        raise OSError(f"Failed to prepare output directory: {e}")

    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        try:
            pbar.set_description(f"Processing {tre_ID}")
            t = Tree(tre_path)
            t1 = scale_support(t, scale)
            tree_str = t1.write(format=0)
            write_tree_to_newick(tree_str, tre_ID, dir_path)
        except Exception as e:
            print(f"Error processing {tre_ID}: {e}")
        pbar.update(1)
    pbar.close()

if __name__ == "__main__":
    """
    Main entry for scaling support values of trees in batch mode.
    Reads a tree dictionary from a file and writes scaled Newick trees to disk.
    """
    try:
        scale = 100  # 可根据需要修改或通过命令行参数传入
        input_file = '100_nosingle_GF_list.txt'
        tre_dic = read_and_return_dict(input_file)
        support_scaler_main(tre_dic, scale)
    except Exception as e:
        print(f"Error in main execution: {e}")
