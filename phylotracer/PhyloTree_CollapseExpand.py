"""
Phylogenetic tree collapse and optional polytomy expansion for the PhyloTracer pipeline.

This module provides utilities to simplify phylogenetic trees by collapsing
low-support internal nodes and to optionally resolve polytomies after collapsing.
"""

import logging
import os
import shutil

logger = logging.getLogger(__name__)

from tqdm import tqdm

from phylotracer import read_and_return_dict, read_tree
from phylotracer.BranchLength_NumericConverter import write_tree_to_newick

# =========================
# Low-Level Tree Utilities
# =========================


def collapse_nodes(phylo_tree: object, node_support_value: int) -> object:
    """Collapse low-support internal nodes in a phylogenetic tree.

    This function removes internal nodes whose support value is below a specified
    threshold while preserving branch lengths, which is commonly used to simplify
    poorly supported topologies before downstream analysis.

    Args:
        phylo_tree (object): Phylogenetic tree object that provides ``copy`` and
            ``traverse`` methods and includes ``support`` on internal nodes.
        node_support_value (int): Minimum support threshold required to retain
            an internal node.

    Returns:
        object: A modified copy of the input tree with low-support nodes collapsed.

    Raises:
        AttributeError: If the input object lacks required tree methods.

    Assumptions:
        Support values are stored as ``node.support`` and represent confidence
        in the corresponding internal node (e.g., bootstrap support).
    """
    if not hasattr(phylo_tree, "copy") or not hasattr(phylo_tree, "traverse"):
        raise AttributeError("phylo_tree must have 'copy' and 'traverse' methods.")
    phylo_tree_copy = phylo_tree.copy()
    for node in phylo_tree_copy.traverse():
        if not node.is_leaf() and not node.is_root():
            if getattr(node, "support", 0) < node_support_value:
                node.delete(preserve_branch_length=True)
    return phylo_tree_copy


# =========================
# Batch Processing & Orchestration
# =========================


def collapse_expand_main(
    tree_dict: dict,
    support_value: int,
    revert: bool = False,
) -> None:
    """Batch-collapse low-support nodes across multiple trees.

    This routine reads multiple trees, collapses internal nodes below a support
    threshold, and optionally resolves polytomies to produce a deterministic
    bifurcating topology for visualization and downstream analyses.

    Args:
        tree_dict (dict): Mapping from tree identifiers to Newick file paths.
        support_value (int): Minimum support threshold for retaining internal nodes.
        revert (bool): If True, resolves polytomies and sorts descendants after
            collapsing nodes.

    Raises:
        ValueError: If ``tree_dict`` is not a dictionary.
        OSError: If output directory creation or cleanup fails.

    Assumptions:
        Tree paths point to valid Newick files and identifiers are unique.
    """
    if not isinstance(tree_dict, dict):
        raise ValueError("tree_dict must be a dictionary.")

    cwd = os.getcwd()
    default_dir = "collapse_expand_tree"
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
            collapsed_tree = collapse_nodes(tree, support_value)
            if revert:
                collapsed_tree.resolve_polytomy(recursive=True)
                collapsed_tree.sort_descendants("support")
            tree_str = collapsed_tree.write(format=0)
            write_tree_to_newick(tree_str, tree_id, dir_path)
        except Exception as exc:
            logger.error("Error processing %s: %s", tree_id, exc)
        pbar.update(1)
    pbar.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Collapse or expand low-support nodes in trees.",
    )
    parser.add_argument(
        "input_file",
        help="Path to the tree list file.",
    )
    parser.add_argument(
        "-s", "--support-value",
        type=int,
        default=50,
        help="Minimum support threshold for retaining internal nodes (default: 50).",
    )
    parser.add_argument(
        "--revert",
        action="store_true",
        default=False,
        help="If set, resolve polytomies instead of collapsing nodes.",
    )
    args = parser.parse_args()

    tre_dic = read_and_return_dict(args.input_file)
    collapse_expand_main(tre_dic, args.support_value, args.revert)
