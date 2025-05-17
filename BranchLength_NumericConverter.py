from __init__ import *

def trans_branch_length(Phylo_t: object, decimal_places: int = 10) -> str:
    """
    Convert a phylogenetic tree object to a Newick string with branch lengths formatted to a specified number of decimal places.

    Args:
        Phylo_t (object): The phylogenetic tree object, must have a 'write' method.
        decimal_places (int): Number of decimal places for branch lengths (default is 10).

    Returns:
        str: Newick string representation of the tree with formatted branch lengths.

    Raises:
        ValueError: If decimal_places is not a non-negative integer.
        AttributeError: If Phylo_t does not have a 'write' method.
    """
    if not isinstance(decimal_places, int) or decimal_places < 0:
        raise ValueError("decimal_places must be a non-negative integer.")
    if not hasattr(Phylo_t, "write"):
        raise AttributeError("Phylo_t must have a 'write' method.")
    try:
        tree_str = Phylo_t.write(format=0, dist_formatter='%.{}f'.format(decimal_places))
    except Exception as e:
        raise RuntimeError(f"Failed to convert tree to Newick string: {e}")
    return tree_str

def write_tree_to_newick(tree_str: str, tre_ID: str, dir_path: str) -> None:
    """
    Write a Newick tree string to a file in the specified directory.

    Args:
        tree_str (str): The Newick tree string.
        tre_ID (str): The tree identifier, used as the file name.
        dir_path (str): The directory path to save the file.

    Raises:
        ValueError: If tree_str is not a string.
        OSError: If file writing fails.
    """
    if not isinstance(tree_str, str):
        raise ValueError("tree_str must be a string.")
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    try:
        with open(os.path.join(dir_path, tre_ID + '.nwk'), 'w') as f:
            f.write(tree_str + '\n')
    except Exception as e:
        raise OSError(f"Failed to write tree to file: {e}")

def branch_length_numeric_converter_main(tre_dic: dict, decimal_places: int = 10) -> None:
    """
    Convert branch lengths of multiple trees to a specified decimal precision and write to files.

    Args:
        tre_dic (dict): Dictionary mapping tree IDs to tree file paths.
        decimal_places (int): Number of decimal places for branch lengths.

    Raises:
        ValueError: If tre_dic is not a dict or decimal_places is invalid.
        OSError: If directory or file operations fail.
    """
    if not isinstance(tre_dic, dict):
        raise ValueError("tre_dic must be a dictionary.")
    if decimal_places is not None and (not isinstance(decimal_places, int) or decimal_places < 0):
        raise ValueError("decimal_places must be a non-negative integer or None.")

    dir_path = os.path.join(os.getcwd(), "converter_tree/")
    try:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path)
    except Exception as e:
        raise OSError(f"Failed to prepare output directory: {e}")

    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        try:
            t = Tree(tre_path)
            if decimal_places is None:
                tree_str = trans_branch_length(t)
            else:
                tree_str = trans_branch_length(t, decimal_places)
            write_tree_to_newick(tree_str, tre_ID, dir_path)
        except Exception as e:
            print(f"Error processing {tre_ID}: {e}")
        pbar.update(1)
    pbar.close()

if __name__ == "__main__":
    """
    Main entry for converting branch lengths of trees in batch mode.
    Reads a tree dictionary from a file and writes formatted Newick trees to disk.
    """
    try:
        input_file = 'test.txt'
        decimal_places = 10  # 可以根据需要修改或通过命令行参数传入
        tre_dic = read_and_return_dict(input_file)
        branch_length_numeric_converter_main(tre_dic, decimal_places)
    except Exception as e:
        print(f"Error in main execution: {e}")
    



