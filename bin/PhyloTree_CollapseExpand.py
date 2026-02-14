"""
Phylogenetic tree collapse and optional polytomy expansion for the PhyloTracer pipeline.

This module provides utilities to simplify phylogenetic trees by collapsing
low-support internal nodes and to optionally resolve polytomies after collapsing.
"""

from __init__ import *
from BranchLength_NumericConverter import write_tree_to_newick

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
        if not node.is_leaf():
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

    dir_path = os.path.join(os.getcwd(), "collapse_expand_tree/")
    try:
        if os.path.exists(dir_path):
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
            print(f"Error processing {tree_id}: {exc}")
        pbar.update(1)
    pbar.close()


if __name__ == "__main__":
    """
    Main entry for collapsing nodes in trees in batch mode.

    This execution path reads a tree dictionary from a file and writes the
    collapsed Newick trees to disk.
    """
    try:
        tre_dic = read_and_return_dict("GF.txt")
        support_value = 50
        revert = False  # Modify as needed.
        collapse_expand_main(tre_dic, support_value, revert)
    except Exception as e:
        print(f"Error in main execution: {e}")
