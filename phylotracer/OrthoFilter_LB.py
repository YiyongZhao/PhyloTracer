"""
Long-branch pruning and visualization for gene trees in the PhyloTracer pipeline.

This module identifies outlier terminal branches in gene trees, optionally
produces before/after visualizations, and writes pruned trees for downstream
orthology-aware analyses.
"""

import os
import shutil
from typing import Dict

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from ete3 import NodeStyle, TextFace, Tree, TreeStyle
from pypdf import PdfReader, PdfWriter
from tqdm import tqdm

from phylotracer import (
    gene_id_transfer,
    get_species_set,
    num_tre_node,
    read_and_return_dict,
    rename_input_tre,
)
from phylotracer.BranchLength_NumericConverter import (
    trans_branch_length,
    write_tree_to_newick,
)

# ======================================================
# Section 1: Tree Property Utilities
# ======================================================


def has_multiple_copies(tree: object) -> bool:
    """
    Determine whether a gene tree contains multiple gene copies.

    Parameters
    ----------
    tree : object
        ETE tree object representing a gene tree.

    Returns
    -------
    bool
        True if the number of leaves differs from the number of unique species.

    Assumptions
    -----------
    Leaf names encode species identifiers that can be parsed by
    ``get_species_set`` from the project utilities.
    """
    return len(tree.get_leaf_names()) != len(get_species_set(tree))


def get_average_tip_length(tree: object) -> float:
    """
    Compute the mean branch length of terminal leaves.

    Parameters
    ----------
    tree : object
        ETE tree object representing a gene tree.

    Returns
    -------
    float
        Mean terminal branch length across all leaves.

    Assumptions
    -----------
    The tree contains at least one leaf and branch lengths are defined.
    """
    return sum(leaf.dist for leaf in tree) / len(tree)


def get_average_node_length(subtree: object) -> float:
    """
    Compute the mean distance from a subtree root to its descendant leaves.

    Parameters
    ----------
    subtree : object
        ETE tree object representing a subtree.

    Returns
    -------
    float
        Mean distance from the subtree root to each descendant leaf.

    Assumptions
    -----------
    Branch lengths are defined for all relevant edges.
    """
    total_distance = sum(
        subtree.get_distance(leaf) + subtree.dist for leaf in subtree
    )
    return total_distance / len(subtree)


# ======================================================
# Section 2: Visualization Utilities
# ======================================================

def create_color_mapping(voucher2taxa: Dict[str, str]) -> Dict[str, str]:
    """
    Assign distinct colors to taxa labels for tree visualization.

    Parameters
    ----------
    voucher2taxa : Dict[str, str]
        Mapping from voucher identifiers to taxa names.

    Returns
    -------
    Dict[str, str]
        Mapping from voucher identifiers to color-annotated taxa labels.

    Assumptions
    -----------
    Taxa names are finite and can be uniquely mapped to colors.
    """
    unique_taxa = sorted(set(voucher2taxa.values()))
    cmap = plt.get_cmap("gist_rainbow")
    color_values = [
        mcolors.rgb2hex(cmap(i))
        for i in np.linspace(0, 1, len(unique_taxa))
    ]
    taxa_color = dict(zip(unique_taxa, color_values))
    return {k: f"{v}*{taxa_color[v]}" for k, v in voucher2taxa.items()}


def style_tree(
    tree: object,
    color_map: Dict[str, str],
    renamed2gene: Dict[str, str],
) -> object:
    """
    Apply node styles and colored gene labels to a tree.

    Parameters
    ----------
    tree : object
        ETE tree object to be styled.
    color_map : Dict[str, str]
        Mapping from species identifiers to color-annotated taxa labels.
    renamed2gene : Dict[str, str]
        Mapping from renamed gene identifiers to original gene names.

    Returns
    -------
    object
        Styled ETE tree object.

    Assumptions
    -----------
    Leaf names encode species identifiers separated by underscores.
    """
    for node in tree.traverse():
        style = NodeStyle()
        style["size"] = 0
        style["shape"] = "circle"
        style["fgcolor"] = "black"
        node.set_style(style)

        if node.is_leaf():
            species = node.name.split("_")[0]
            gene_name = renamed2gene.get(node.name, node.name)

            if species in color_map:
                color = color_map[species].split("*")[-1]
                node.add_face(
                    TextFace(gene_name, fgcolor=color, fstyle="italic"),
                    column=0,
                )
    return tree


def generate_tree_pdf(tree_id: str, tree: object, suffix: str) -> None:
    """
    Render a tree into a single-page PDF file.

    Parameters
    ----------
    tree_id : str
        Identifier of the gene tree.
    tree : object
        ETE tree object to be rendered.
    suffix : str
        Suffix appended to the output filename for disambiguation.

    Returns
    -------
    None

    Assumptions
    -----------
    The ETE rendering backend is available in the runtime environment.
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(f"{tree_id}_{suffix}", fsize=10), column=0)
    tree.render(file_name=f"{tree_id}_{suffix}.pdf", tree_style=ts)


def merge_pdfs(left_pdf: str, right_pdf: str, output_pdf: str) -> None:
    """
    Merge two single-page PDFs into a single horizontal layout.

    Parameters
    ----------
    left_pdf : str
        Path to the left PDF file.
    right_pdf : str
        Path to the right PDF file.
    output_pdf : str
        Path to the merged output PDF.

    Returns
    -------
    None

    Assumptions
    -----------
    Each input PDF contains at least one page.
    """
    writer = PdfWriter()

    with open(left_pdf, "rb") as f1, open(right_pdf, "rb") as f2:
        pdf1, pdf2 = PdfReader(f1), PdfReader(f2)
        page1, page2 = pdf1.pages[0], pdf2.pages[0]

        width = float(page1.mediabox.width) + float(page2.mediabox.width)
        height = max(float(page1.mediabox.height), float(page2.mediabox.height))

        new_page = writer.add_blank_page(width, height)
        new_page.merge_translated_page(page1, 0, 0)
        new_page.merge_translated_page(page2, float(page1.mediabox.width), 0)

        with open(output_pdf, "wb") as out:
            writer.write(out)


# ======================================================
# Section 3: Long Branch Pruning Logic (Core)
# ======================================================

def remove_long_branches(
    tree: object,
    abs_threshold: float,
    rel_threshold: float,
    log_handle,
    tree_id: str,
    renamed2gene: Dict[str, str],
) -> object:
    """
    Prune leaves whose branch lengths exceed absolute and relative thresholds.

    Parameters
    ----------
    tree : object
        ETE tree object representing a gene tree.
    abs_threshold : float
        Absolute threshold on the root-relative branch ratio.
    rel_threshold : float
        Relative threshold on the sister-relative branch ratio.
    log_handle : object
        File handle for recording pruning decisions.
    tree_id : str
        Identifier of the gene tree.
    renamed2gene : Dict[str, str]
        Mapping from renamed gene identifiers to original gene names.

    Returns
    -------
    object
        Pruned ETE tree object.

    Assumptions
    -----------
    Tree is rooted and branch lengths are comparable across leaves.
    """
    pruned_tree = tree.copy()
    to_remove: Set[str] = set()

    avg_tip_length = get_average_tip_length(pruned_tree)

    for leaf in pruned_tree:
        leaf_name = leaf.name
        leaf_dist = leaf.dist
        gene_name = renamed2gene.get(leaf_name, leaf_name)

        if leaf_dist == 0:
            log_handle.write(f"{tree_id}\t\t{gene_name}\t0\t0\n")
            continue

        root_ratio = (leaf_dist - avg_tip_length) / avg_tip_length

        sisters = leaf.get_sisters()
        if not sisters:
            log_handle.write(f"{tree_id}\t\t{gene_name}\t0\t0\n")
            continue
        sister = sisters[0]
        if sister.is_leaf():
            sister_avg = sister.dist or 1e-6
        else:
            sister_avg = get_average_node_length(sister) or 1e-6

        sister_ratio = (leaf_dist - sister_avg) / sister_avg

        if root_ratio >= abs_threshold:
            if sister_ratio >= rel_threshold:
                log_handle.write(
                    f"{tree_id}\t*\t{gene_name}\t{root_ratio}\t{sister_ratio}\n"
                )
                to_remove.add(leaf_name)
            else:
                log_handle.write(
                    f"{tree_id}\t\t{gene_name}\t{root_ratio}\t{sister_ratio}\n"
                )
        else:
            log_handle.write(
                f"{tree_id}\t\t{gene_name}\t{root_ratio}\t{sister_ratio}\n"
            )

    keep_leaves = set(pruned_tree.get_leaf_names()) - to_remove
    pruned_tree.prune(keep_leaves, preserve_branch_length=True)
    return pruned_tree


# ======================================================
# Section 4: Main Pipeline (Orchestration)
# ======================================================

def prune_main_LB(
    tree_dict: Dict[str, str],
    voucher2taxa: Dict[str, str],
    gene2renamed: Dict[str, str],
    renamed2gene: Dict[str, str],
    absolute_branch_length: float = 5,
    relative_branch_length: float = 5,
    visual: bool = False,
) -> None:
    """
    Run the long-branch pruning workflow across a set of gene trees.

    Parameters
    ----------
    tree_dict : Dict[str, str]
        Mapping from tree identifiers to Newick file paths.
    voucher2taxa : Dict[str, str]
        Mapping from voucher identifiers to taxa names.
    gene2renamed : Dict[str, str]
        Mapping from original gene identifiers to renamed identifiers.
    renamed2gene : Dict[str, str]
        Mapping from renamed gene identifiers to original gene names.
    absolute_branch_length : float, optional
        Absolute threshold for the root-relative branch ratio.
    relative_branch_length : float, optional
        Relative threshold for the sister-relative branch ratio.
    visual : bool, optional
        Whether to generate before/after tree visualization PDFs.

    Returns
    -------
    None

    Assumptions
    -----------
    Input trees are valid Newick files and names can be mapped by the provided
    gene identifier dictionaries.
    """
    color_map = create_color_mapping(voucher2taxa)

    base_dir = os.getcwd()
    pruned_dir = os.path.join(base_dir, "orthofilter_lb/pruned_tree")
    log_dir = os.path.join(base_dir, "orthofilter_lb/long_branch_gene")
    if visual:
        pdf_dir = os.path.join(base_dir, "orthofilter_lb/pruned_tree_pdf")
    else:
        pdf_dir = None


    for d in (pruned_dir, log_dir, pdf_dir):
        if d is not None:
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d, exist_ok=True)

    pbar = tqdm(tree_dict.items(), desc="Processing trees", unit="tree")

    for tree_id, tree_path in pbar:
        pbar.set_description(f"Processing {tree_id}")
        tree = Tree(tree_path)
        tree.ladderize()
        tree.resolve_polytomy(recursive=True)
        tree.sort_descendants("support")

        tree = rename_input_tre(tree, gene2renamed)
        num_tre_node(tree)

        log_path = os.path.join(log_dir, f"{tree_id}_delete_gene.txt")
        with open(log_path, "w") as log:
            log.write(
                "tre_ID\tlong_branch_label\tgene\t"
                "root_relative_branch_ratio\tsister_relative_branch_ratio\n"
            )

            if visual:
                styled = style_tree(tree, color_map, renamed2gene)
                generate_tree_pdf(tree_id, styled, "before")

            pruned_tree = remove_long_branches(
                tree,
                absolute_branch_length,
                relative_branch_length,
                log,
                tree_id,
                renamed2gene,
            )

        if visual:
            generate_tree_pdf(tree_id, pruned_tree, "after")
            merge_pdfs(
                f"{tree_id}_before.pdf",
                f"{tree_id}_after.pdf",
                os.path.join(pdf_dir, f"{tree_id}.pdf"),
            )
            os.remove(f"{tree_id}_before.pdf")
            os.remove(f"{tree_id}_after.pdf")

        restored_tree = rename_input_tre(pruned_tree, renamed2gene)
        tree_str = trans_branch_length(restored_tree)
        write_tree_to_newick(tree_str, tree_id, pruned_dir)
    pbar.close()


# ======================================================
# Section 5: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Long-branch pruning filter")
    parser.add_argument("--input_GF_list", required=True, help="Gene family list file")
    parser.add_argument("--input_imap", required=True, help="Imap file")
    parser.add_argument("--absolute_branch_length", type=float, default=5, help="Absolute branch length cutoff")
    parser.add_argument("--relative_branch_length", type=float, default=5, help="Relative branch length cutoff")
    parser.add_argument("--visual", action="store_true", help="Enable visual output")
    args = parser.parse_args()

    tree_dict = read_and_return_dict(args.input_GF_list)
    (
        gene2renamed,
        renamed2gene,
        voucher2taxa,
        taxa2voucher,
    ) = gene_id_transfer(args.input_imap)

    prune_main_LB(
        tree_dict,
        voucher2taxa,
        gene2renamed,
        renamed2gene,
        absolute_branch_length=args.absolute_branch_length,
        relative_branch_length=args.relative_branch_length,
        visual=args.visual,
    )
