"""
Support-value scaling utilities for the PhyloTracer pipeline.

This module rescales node support values to a consistent numeric range and
writes standardized Newick outputs for downstream analyses.
"""

import logging
import os
import shutil

logger = logging.getLogger(__name__)

from ete3 import Tree
from tqdm import tqdm

from phylotracer import read_and_return_dict, read_tree
from phylotracer.BranchLength_NumericConverter import write_tree_to_newick

# =========================
# Low-Level Tree Utilities
# =========================


def scale_support(phylo_tree: object, scale: str = "100") -> object:
    """Scale support values across all nodes in a phylogenetic tree.

    Args:
        phylo_tree (object): Phylogenetic tree object with ``traverse`` and
            ``support`` attributes.
        scale (str or int): Target scale; ``100`` scales values to [1, 100],
            otherwise values are scaled down to [0, 1] when appropriate.

    Returns:
        object: The same tree object with rescaled support values.

    Raises:
        AttributeError: If the tree does not provide a ``traverse`` method.
        RuntimeError: If support values cannot be evaluated.

    Assumptions:
        Support values are either in [0, 1] or [0, 100], and rescaling preserves
        their relative magnitudes.
    """
    if not hasattr(phylo_tree, "traverse"):
        raise AttributeError("phylo_tree must have a 'traverse' method.")
    try:
        max_support = max(getattr(node, "support", 0) for node in phylo_tree.traverse())
    except Exception as exc:
        raise RuntimeError(f"Failed to get support values: {exc}")

    if isinstance(scale, str) and scale.isdigit():
        scale = int(scale)

    if max_support < 1:
        if scale == 100:
            for node in phylo_tree.traverse():
                if hasattr(node, "support"):
                    node.support *= 100
        else:
            logger.info(
                "No need to scale down support, as all tree node supports are already < 1."
            )
    else:
        if scale != 100:
            for node in phylo_tree.traverse():
                if hasattr(node, "support"):
                    node.support /= 100
        else:
            logger.info(
                "No need to scale support, as tree node supports are already in the range [1, 100]."
            )
    return phylo_tree


# =========================
# Batch Processing & Orchestration
# =========================


def support_scaler_main(tree_dict: dict, scale: str = "100") -> None:
    """Batch-scale support values and write Newick outputs.

    Args:
        tree_dict (dict): Mapping from tree identifiers to tree file paths.
        scale (str or int): Target scale for support values.

    Raises:
        ValueError: If ``tree_dict`` is not a dictionary.
        OSError: If output directory creation or cleanup fails.

    Assumptions:
        Tree paths are valid and refer to Newick-readable trees.
    """
    if not isinstance(tree_dict, dict):
        raise ValueError("tree_dict must be a dictionary.")

    cwd = os.getcwd()
    default_dir = "support_scaler_tree"
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
            scaled_tree = scale_support(tree, scale)
            tree_str = scaled_tree.write(format=0)
            write_tree_to_newick(tree_str, tree_id, dir_path)
        except Exception as exc:
            logger.error("Error processing %s: %s", tree_id, exc)
        pbar.update(1)
    pbar.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Scale support values of trees in batch mode.",
    )
    parser.add_argument(
        "input_file",
        help="Path to the tree list file.",
    )
    parser.add_argument(
        "-s", "--scale",
        type=int,
        default=100,
        help="Scale factor for support values (default: 100).",
    )
    args = parser.parse_args()

    tre_dic = read_and_return_dict(args.input_file)
    support_scaler_main(tre_dic, args.scale)
