"""
Tree topology summarization and visualization for the PhyloTracer pipeline.

This module groups gene trees by topology, writes absolute and relative
summaries, and optionally renders representative topologies for reporting.
"""
from __future__ import annotations

import logging
import math
import os
import shutil
import tempfile
try:
    import matplotlib as mpl
    mpl.rcParams["pdf.fonttype"] = 42
except Exception:
    mpl = None

logger = logging.getLogger(__name__)

from ete3 import Tree
try:
    from ete3 import NodeStyle, TextFace, TreeStyle
except ImportError:
    NodeStyle = None
    TextFace = None
    TreeStyle = None
try:
    import fitz  # PyMuPDF
except ImportError:
    fitz = None
try:
    from pypdf import PdfReader, PdfWriter
    from pypdf import Transformation
except ImportError:
    PdfReader = None
    PdfWriter = None
    Transformation = None

from phylotracer import (
    gene_id_transfer,
    get_species_set,
    read_and_return_dict,
    read_tree,
    rename_input_tre,
)

# ======================================================
# Section 1: Tree Simplification and Grouping Utilities
# ======================================================


def get_only_sps_tree(Phylo_t: object) -> object:
    """
    Simplify leaf names to species identifiers only.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be simplified.

    Returns
    -------
    object
        A tree with leaf names reduced to species codes.

    Assumptions
    -----------
    Leaf names encode species identifiers before the first underscore.
    """
    Phylo_t_c = Phylo_t.copy()
    for node in Phylo_t_c:
        node.name = node.name.split("_")[0]
    return Phylo_t_c


def get_max_tree(trees: list) -> object:
    """
    Identify the tree with the maximum number of leaves.

    Parameters
    ----------
    trees : list[object]
        List of ETE tree objects.

    Returns
    -------
    object
        Tree with the highest leaf count, or None if input is empty.

    Assumptions
    -----------
    Leaf counts are computed consistently by ``get_leaves``.
    """
    max_leaf_count = 0
    max_leaf_count_tree = None

    for tree in trees:
        leaf_count = len(tree.get_leaves())
        if leaf_count > max_leaf_count:
            max_leaf_count = leaf_count
            max_leaf_count_tree = tree

    return max_leaf_count_tree


def group_trees_by_topology_with_ids(
    trees_with_ids: list[tuple],
    result_dict: dict,
) -> list[tuple]:
    """
    Group trees by topology using Robinson-Foulds distance.

    Parameters
    ----------
    trees_with_ids : list[tuple]
        List of (tree_id, tree_object) tuples.
    result_dict : dict
        Dictionary that stores grouped trees by Newick string keys.

    Returns
    -------
    list[tuple]
        Trees that do not match the topology of the first tree in the list.

    Assumptions
    -----------
    RF distance of zero indicates identical topology.
    """
    if len(trees_with_ids) == 1:
        tree_id, tree = trees_with_ids[0]
        result_dict[tree.write(format=9)] = [trees_with_ids[0]]
        return []
    else:
        same_rf_trees = []
        different_trees = []
        first_tree_id, first_tree = trees_with_ids[0]

        for tree_id, tree in trees_with_ids:
            rf = first_tree.robinson_foulds(tree)[0]
            if rf == 0:
                same_rf_trees.append((tree_id, tree))
            else:
                different_trees.append((tree_id, tree))

        max_tree_tuple = max(same_rf_trees, key=lambda x: len(x[1].get_leaves()))
        result_dict[max_tree_tuple[1].write(format=9)] = same_rf_trees
        return different_trees


def process_tree_with_ids(trees_with_ids: list[tuple], result_dict: dict) -> None:
    """
    Iteratively group trees by topology until all are assigned.

    Parameters
    ----------
    trees_with_ids : list[tuple]
        List of (tree_id, tree_object) tuples.
    result_dict : dict
        Dictionary that stores grouped trees by Newick string keys.

    Returns
    -------
    None

    Assumptions
    -----------
    The grouping function partitions the list without loss.
    """
    remaining_trees = trees_with_ids
    while len(remaining_trees) >= 1:
        remaining_trees = group_trees_by_topology_with_ids(
            remaining_trees,
            result_dict,
        )


# ======================================================
# Section 2: Summary Writing Utilities
# ======================================================


def write_relative_summary_with_ids(
    outfile: str,
    dic: dict,
    voucher2taxa_dic: dict,
    top_n: int = None,
) -> None:
    """
    Write relative topology summaries and optionally visualize top ranks.

    Parameters
    ----------
    outfile : str
        Base name for output files.
    dic : dict
        Mapping from topology Newick strings to lists of (tree_id, tree) tuples.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels for renaming.
    top_n : int, optional
        Number of top topologies to visualize.

    Returns
    -------
    None

    Assumptions
    -----------
    Topology frequency is determined by list length per entry.
    """
    sorted_dict = dict(sorted(dic.items(), key=lambda item: len(item[1]), reverse=True))

    if top_n:
        output_pdf = f"merge_relative_top{top_n}.pdf"
        visualize_top_trees(sorted_dict, output_pdf, voucher2taxa_dic, top_n)

    with open(f"relative_{outfile}.txt", "w") as file:
        file.write("Topology_ID\tTopology_num\tTopology\tTree_IDs\n")
        for index, clade in enumerate(sorted_dict.keys()):
            trees_with_ids = sorted_dict[clade]

            representative_tree_tuple = max(
                trees_with_ids,
                key=lambda x: len(x[1].get_leaves()),
            )
            _, representative_tree = representative_tree_tuple

            tree_renamed = rename_input_tre(representative_tree, voucher2taxa_dic)
            tree_str = tree_renamed.write(format=9)

            tree_ids = [tree_id for tree_id, _ in trees_with_ids]
            tree_ids_str = ",".join(tree_ids)

            file.write(
                f"{index}\t{len(trees_with_ids)}\t{tree_str}\t{tree_ids_str}\n"
            )


def write_absolute_summary(
    outfile: str,
    trees_with_ids: list[tuple],
    voucher2taxa_dic: dict,
    top_n: int = None,
) -> None:
    """
    Summarize absolute topologies and optionally visualize top ranks.

    Parameters
    ----------
    outfile : str
        Output file base name.
    trees_with_ids : list[tuple]
        List of (tree_id, tree_object) tuples.
    voucher2taxa_dic : dict
        Mapping for renaming tree leaves.
    top_n : int, optional
        Number of top topologies to visualize.

    Returns
    -------
    None

    Assumptions
    -----------
    Absolute topology is defined by Newick string identity.
    """
    topolo_dic = {}
    for tree_id, tree in trees_with_ids:
        tree_str = tree.write(format=9)
        if tree_str in topolo_dic:
            topolo_dic[tree_str].append((tree_id, tree))
        else:
            topolo_dic[tree_str] = [(tree_id, tree)]

    sorted_dict = dict(
        sorted(
            topolo_dic.items(),
            key=lambda item: (len(item[1]), len(item[0])),
            reverse=True,
        )
    )

    if top_n:
        output_pdf = f"merge_absolutely_top{top_n}.pdf"
        visualize_top_trees(sorted_dict, output_pdf, voucher2taxa_dic, top_n)

    with open(f"absolute_{outfile}.txt", "w") as file:
        file.write("Topology_ID\tTopology_num\tTopology\tTree_IDs\n")
        for num, k in enumerate(sorted_dict):
            t1 = Tree(k)
            t2 = rename_input_tre(t1, voucher2taxa_dic)
            t1_str = t2.write(format=9)
            tree_ids = [tree_id for tree_id, _ in sorted_dict[k]]
            tree_ids_str = ",".join(tree_ids)
            file.write(
                f"{num}\t{len(sorted_dict[k])}\t{t1_str}\t{tree_ids_str}\n"
            )


# ======================================================
# Section 3: Topology Visualization
# ======================================================


def visualize_top_trees(
    tree_count_dict: dict,
    output_path: str,
    voucher2taxa_dic: dict,
    top_n: int = 10,
) -> None:
    """
    Render the top N topologies as a combined grid image.

    Parameters
    ----------
    tree_count_dict : dict
        Mapping from topology Newick strings to lists of trees.
    output_path : str
        Output path for the combined image.
    voucher2taxa_dic : dict
        Mapping for renaming tree leaves.
    top_n : int, optional
        Number of top topologies to render.

    Returns
    -------
    None

    Assumptions
    -----------
    Tree rendering is available, and PDF merge backend (PyMuPDF or pypdf)
    is installed.
    """
    def _layout_fixed_leaf_font(node):
        if node.is_leaf():
            # Nature-like panel typography baseline.
            node.add_face(TextFace(node.name, fsize=6), column=0, position="branch-right")

    def _get_pdf_panel_sizes(pdf_files: list[str]) -> tuple[float, float]:
        max_w = 0.0
        max_h = 0.0
        if fitz is not None:
            for path in pdf_files:
                src = fitz.open(path)
                rect = src[0].rect
                max_w = max(max_w, float(rect.width))
                max_h = max(max_h, float(rect.height))
                src.close()
            return max_w, max_h
        if PdfReader is not None:
            for path in pdf_files:
                with open(path, "rb") as handle:
                    p = PdfReader(handle).pages[0]
                    max_w = max(max_w, float(p.mediabox.width))
                    max_h = max(max_h, float(p.mediabox.height))
            return max_w, max_h
        return 0.0, 0.0

    def _merge_pdf_with_pymupdf(
        pdf_files: list[str],
        output_pdf_path: str,
        cols: int,
        rows: int,
        cell_w: float,
        cell_h: float,
        gap_x: float,
        gap_y: float,
        margin: float,
    ) -> bool:
        if fitz is None:
            return False
        # Keep panel size unchanged to avoid any compression/scaling.
        page_width = 2 * margin + cols * cell_w + (cols - 1) * gap_x
        page_height = 2 * margin + rows * cell_h + (rows - 1) * gap_y
        merged_doc = fitz.open()
        page = merged_doc.new_page(width=page_width, height=page_height)
        x = margin
        y = margin
        count = 0
        for panel_pdf in pdf_files:
            src = fitz.open(panel_pdf)
            src_rect = src[0].rect
            # Zero-scaling placement: keep original panel width/height.
            rect = fitz.Rect(x, y, x + float(src_rect.width), y + float(src_rect.height))
            page.show_pdf_page(rect, src, 0)
            src.close()

            count += 1
            if count % cols == 0:
                x = margin
                y += cell_h + gap_y
            else:
                x += cell_w + gap_x

            if count == rows * cols and panel_pdf != pdf_files[-1]:
                page = merged_doc.new_page(width=page_width, height=page_height)
                x = margin
                y = margin
                count = 0

        merged_doc.save(output_pdf_path)
        merged_doc.close()
        return True

    def _merge_pdf_with_pypdf(
        pdf_files: list[str],
        output_pdf_path: str,
        cols: int,
        rows: int,
        cell_w: float,
        cell_h: float,
        gap_x: float,
        gap_y: float,
        margin: float,
    ) -> bool:
        if PdfReader is None or PdfWriter is None or Transformation is None:
            return False
        page_width = 2 * margin + cols * cell_w + (cols - 1) * gap_x
        page_height = 2 * margin + rows * cell_h + (rows - 1) * gap_y
        writer = PdfWriter()
        page = writer.add_blank_page(width=page_width, height=page_height)
        x = margin
        y = margin
        count = 0
        for panel_pdf in pdf_files:
            with open(panel_pdf, "rb") as handle:
                src_page = PdfReader(handle).pages[0]
                # Zero-scaling placement for publication fidelity.
                page.merge_translated_page(src_page, tx=x, ty=y)
            count += 1
            if count % cols == 0:
                x = margin
                y += cell_h + gap_y
            else:
                x += cell_w + gap_x
            if count == rows * cols and panel_pdf != pdf_files[-1]:
                page = writer.add_blank_page(width=page_width, height=page_height)
                x = margin
                y = margin
                count = 0
        with open(output_pdf_path, "wb") as out_handle:
            writer.write(out_handle)
        return True

    cols = math.ceil(math.sqrt(top_n))
    rows = math.ceil(top_n / cols)
    panel_h = 2200
    gap_x = 60
    gap_y = 60
    margin = 40
    sorted_trees = list(tree_count_dict.items())[:top_n]
    temp_dir = tempfile.mkdtemp(prefix="temp_trees_")
    tree_panels_pdf = []
    for i, (tree_str, count) in enumerate(sorted_trees):
        tree0 = read_tree(tree_str)
        tree = rename_input_tre(tree0, voucher2taxa_dic)
        tree.ladderize()
        tree.resolve_polytomy(recursive=True)
        tree.sort_descendants("support")

        tree.convert_to_ultrametric()

        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 0.5
        nstyle["hz_line_width"] = 0.5
        nstyle["vt_line_type"] = 0
        nstyle["hz_line_type"] = 0
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        for node in tree.traverse():
            node.set_style(nstyle)
        ts = TreeStyle()
        ts.scale = 10
        ts.show_leaf_name = False
        ts.layout_fn = _layout_fixed_leaf_font
        ts.title.add_face(TextFace(f"Count: {len(count)}", fsize=6), column=1)
        ts.show_scale = False
        ts.force_topology = True
        ts.branch_vertical_margin = 4
        ts.margin_left = 20
        ts.margin_right = 20
        ts.margin_top = 20
        ts.margin_bottom = 20
        panel_pdf = os.path.join(temp_dir, f"tree_{i + 1}.pdf")
        # Keep natural width; only fix height to avoid width/height squeeze.
        tree.render(panel_pdf, h=panel_h, tree_style=ts)
        tree_panels_pdf.append(panel_pdf)
    if tree_panels_pdf:
        output_pdf = output_path if output_path.lower().endswith(".pdf") else f"{os.path.splitext(output_path)[0]}.pdf"
        cell_w, cell_h = _get_pdf_panel_sizes(tree_panels_pdf)
        if cell_w <= 0 or cell_h <= 0:
            raise RuntimeError("Unable to read panel PDF size for vector merge.")
        merged_ok = _merge_pdf_with_pymupdf(
            tree_panels_pdf, output_pdf, cols, rows, cell_w, cell_h, gap_x, gap_y, margin
        )
        if not merged_ok:
            merged_ok = _merge_pdf_with_pypdf(
                tree_panels_pdf, output_pdf, cols, rows, cell_w, cell_h, gap_x, gap_y, margin
            )
        if not merged_ok:
            raise RuntimeError(
                "Failed to merge vector PDF panels: install pymupdf or pypdf."
            )
    shutil.rmtree(temp_dir, ignore_errors=True)


# ======================================================
# Section 4: Main Pipeline (Orchestration)
# ======================================================


def statistical_main(
    tre_dic: dict,
    outfile: str,
    gene2new_named_gene_dic: dict,
    voucher2taxa_dic: dict,
    top_n: int,
) -> None:
    """
    Summarize tree topologies and write absolute and relative outputs.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree identifiers to file paths.
    outfile : str
        Base name for output files.
    gene2new_named_gene_dic : dict
        Mapping from gene identifiers to renamed identifiers.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxon labels.
    top_n : int
        Number of top topologies to visualize.

    Returns
    -------
    None

    Assumptions
    -----------
    Input trees are single-copy after species-only simplification.
    """
    trees_with_ids = []
    for k, v in tre_dic.items():
        t = read_tree(v)
        t1 = rename_input_tre(t, gene2new_named_gene_dic)
        t1.sort_descendants()

        t2 = get_only_sps_tree(t1)
        if len(get_species_set(t2)) == 1:
            continue
        if len(get_species_set(t2)) != len(t2):
            logger.error("Gene tree %s is not a single-copy gene tree!", k)
            logger.error(
                "Number of species: %d, Number of genes: %d",
                len(get_species_set(t2)), len(t2)
            )
            logger.error(
                "This program requires single-copy gene trees as input. "
                "Program terminated."
            )
            exit(1)
        trees_with_ids.append((k, t2))

    write_absolute_summary(outfile, trees_with_ids, voucher2taxa_dic, top_n)

    dic_with_ids = {}
    process_tree_with_ids(trees_with_ids, dic_with_ids)
    write_relative_summary_with_ids(outfile, dic_with_ids, voucher2taxa_dic, top_n)


# ======================================================
# Section 5: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Tree topology summarization")
    parser.add_argument("--input_GF_list", required=True, help="Gene family list file")
    parser.add_argument("--input_imap", required=True, help="Imap file")
    parser.add_argument("--output", default="result", help="Output file prefix")
    parser.add_argument("--top_n", type=int, default=10, help="Number of top topologies to report")
    args = parser.parse_args()

    tre_dic = read_and_return_dict(args.input_GF_list)
    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, _ = (
        gene_id_transfer(args.input_imap)
    )
    statistical_main(tre_dic, args.output, gene2new_named_gene_dic, voucher2taxa_dic, top_n=args.top_n)
