
from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

def collapse_nodes(Phylo_t: object, node_support_value: int) -> object:
    """
    Collapse nodes in a phylogenetic tree whose support value is below a threshold.

    Args:
        Phylo_t (object): The phylogenetic tree object, must have 'copy' and 'traverse' methods.
        node_support_value (int): The support value threshold for collapsing nodes.

    Returns:
        object: The collapsed tree object.
    """
    if not hasattr(Phylo_t, "copy") or not hasattr(Phylo_t, "traverse"):
        raise AttributeError("Phylo_t must have 'copy' and 'traverse' methods.")
    Phylo_t_c = Phylo_t.copy()
    for node in Phylo_t_c.traverse():
        if not node.is_leaf():
            if getattr(node, 'support', 0) < node_support_value:
                node.delete(preserve_branch_length=True)
    return Phylo_t_c

def collapse_expand_main(tre_dic: dict, support_value: int, revert: bool = False) -> None:
    """
    Batch process trees to collapse nodes below a support threshold and optionally expand polytomies.

    Args:
        tre_dic (dict): Dictionary mapping tree IDs to tree file paths.
        support_value (int): The support value threshold for collapsing nodes.
        revert (bool): Whether to resolve polytomies and sort descendants after collapsing.

    Raises:
        ValueError: If tre_dic is not a dict.
        OSError: If directory or file operations fail.
    """
    if not isinstance(tre_dic, dict):
        raise ValueError("tre_dic must be a dictionary.")

    dir_path = os.path.join(os.getcwd(), "collapse_expand_tree/")
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
            t1 = collapse_nodes(t, support_value)
            if revert:
                t1.resolve_polytomy(recursive=True)
                t1.sort_descendants("support")
            tree_str = t1.write(format=0)
            write_tree_to_newick(tree_str, tre_ID, dir_path)
        except Exception as e:
            print(f"Error processing {tre_ID}: {e}")
        pbar.update(1)
    pbar.close()

if __name__ == "__main__":
    """
    Main entry for collapsing nodes in trees in batch mode.
    Reads a tree dictionary from a file and writes collapsed Newick trees to disk.
    """
    try:
        tre_dic = read_and_return_dict('GF.txt')
        support_value = 50
        revert = False  # 可根据需要修改
        collapse_expand_main(tre_dic, support_value, revert)
    except Exception as e:
        print(f"Error in main execution: {e}")
    
