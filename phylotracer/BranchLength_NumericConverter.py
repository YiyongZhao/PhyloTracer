"""
Branch-length formatting and export utilities for the PhyloTracer pipeline.

This module converts branch lengths to fixed decimal precision and writes
standardized Newick outputs for downstream analyses and visualization.
"""

import logging
import os
import shutil

logger = logging.getLogger(__name__)

from tqdm import tqdm

from phylotracer import read_tree, read_and_return_dict

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
    os.makedirs(dir_path, exist_ok=True)
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

    cwd = os.getcwd()
    default_dir = "converter_tree"
    dir_path = (
        cwd
        if os.path.basename(os.path.normpath(cwd)) == default_dir
        else os.path.join(cwd, f"{default_dir}/")
    )
    try:
        if dir_path != cwd and os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path, exist_ok=True)
    except Exception as exc:
        raise OSError(f"Failed to prepare output directory: {exc}")

    pbar = tqdm(total=len(tree_dict), desc="Processing trees", unit="tree")
    for tree_id, tree_path in tree_dict.items():
        try:
            pbar.set_description(f"Processing {tree_id}")
            tree = read_tree(tree_path)
            if decimal_places is None:
                tree_str = trans_branch_length(tree)
            else:
                tree_str = trans_branch_length(tree, decimal_places)
            write_tree_to_newick(tree_str, tree_id, dir_path)
        except Exception as exc:
            logger.error("Error processing %s: %s", tree_id, exc)
        pbar.update(1)
    pbar.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert branch lengths of trees to fixed decimal precision.",
    )
    parser.add_argument(
        "input_file",
        help="Path to the tree list file.",
    )
    parser.add_argument(
        "-d", "--decimal-places",
        type=int,
        default=10,
        help="Number of decimal places for branch lengths (default: 10).",
    )
    args = parser.parse_args()

    tre_dic = read_and_return_dict(args.input_file)
    branch_length_numeric_converter_main(tre_dic, args.decimal_places)
