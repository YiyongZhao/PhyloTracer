"""
Branch-length formatting and export utilities for the PhyloTracer pipeline.

This module converts branch lengths to fixed decimal precision and writes
standardized Newick outputs for downstream analyses and visualization.
"""

from __init__ import *

# =========================
# Low-Level Tree Utilities
# =========================


def trans_branch_length(phylo_tree: object, decimal_places: int = 10) -> str:
    """Convert a phylogenetic tree to a Newick string with formatted branch lengths.

    Args:
        phylo_tree (object): Phylogenetic tree object that provides a ``write`` method.
        decimal_places (int): Number of decimal places for branch lengths.

    Returns:
        str: Newick string representation of the tree with formatted branch lengths.

    Raises:
        ValueError: If ``decimal_places`` is not a non-negative integer.
        AttributeError: If ``phylo_tree`` does not provide a ``write`` method.
        RuntimeError: If tree serialization fails.

    Assumptions:
        The tree object supports ``write(format=0, dist_formatter=...)`` as in ete3.
    """
    if not isinstance(decimal_places, int) or decimal_places < 0:
        raise ValueError("decimal_places must be a non-negative integer.")
    if not hasattr(phylo_tree, "write"):
        raise AttributeError("phylo_tree must have a 'write' method.")
    try:
        tree_str = phylo_tree.write(
            format=0,
            dist_formatter="%.{}f".format(decimal_places),
        )
    except Exception as exc:
        raise RuntimeError(f"Failed to convert tree to Newick string: {exc}")
    return tree_str


def write_tree_to_newick(tree_str: str, tree_id: str, dir_path: str) -> None:
    """Write a Newick tree string to disk.

    Args:
        tree_str (str): Newick tree string to write.
        tree_id (str): Tree identifier used as the output filename stem.
        dir_path (str): Output directory path.

    Raises:
        ValueError: If ``tree_str`` is not a string.
        OSError: If file creation or writing fails.

    Assumptions:
        The output directory is writable; a trailing newline is appended.
    """
    if not isinstance(tree_str, str):
        raise ValueError("tree_str must be a string.")
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    try:
        with open(os.path.join(dir_path, tree_id + ".nwk"), "w") as f:
            f.write(tree_str + "\n")
    except Exception as exc:
        raise OSError(f"Failed to write tree to file: {exc}")


# =========================
# Batch Processing & Orchestration
# =========================


def branch_length_numeric_converter_main(
    tree_dict: dict,
    decimal_places: int = 10,
) -> None:
    """Batch-convert branch lengths to fixed precision and write Newick files.

    Args:
        tree_dict (dict): Mapping from tree identifiers to tree file paths.
        decimal_places (int): Number of decimal places for branch lengths.

    Raises:
        ValueError: If ``tree_dict`` is not a dict or ``decimal_places`` is invalid.
        OSError: If output directory creation or cleanup fails.

    Assumptions:
        Tree paths are valid and refer to Newick-readable trees.
    """
    if not isinstance(tree_dict, dict):
        raise ValueError("tree_dict must be a dictionary.")
    if decimal_places is not None and (
        not isinstance(decimal_places, int) or decimal_places < 0
    ):
        raise ValueError("decimal_places must be a non-negative integer or None.")

    dir_path = os.path.join(os.getcwd(), "converter_tree/")
    try:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path)
    except Exception as exc:
        raise OSError(f"Failed to prepare output directory: {exc}")

    pbar = tqdm(total=len(tree_dict), desc="Processing trees", unit="tree")
    for tree_id, tree_path in tree_dict.items():
        try:
            tree = read_tree(tree_path)
            if decimal_places is None:
                tree_str = trans_branch_length(tree)
            else:
                tree_str = trans_branch_length(tree, decimal_places)
            write_tree_to_newick(tree_str, tree_id, dir_path)
        except Exception as exc:
            print(f"Error processing {tree_id}: {exc}")
        pbar.update(1)
    pbar.close()


if __name__ == "__main__":
    """
    Main entry for converting branch lengths of trees in batch mode.
    Reads a tree dictionary from a file and writes formatted Newick trees to disk.
    """
    try:
        input_file = "test.txt"
        decimal_places = 10  # You can modify it as needed or pass it via command-line parameters.
        tre_dic = read_and_return_dict(input_file)
        branch_length_numeric_converter_main(tre_dic, decimal_places)
    except Exception as e:
        print(f"Error in main execution: {e}")
