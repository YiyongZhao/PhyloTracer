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
    from pypdf import PdfReader, PdfWriter, Transformation
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

logger = logging.getLogger(__name__)

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
            try:
                rf = first_tree.robinson_foulds(tree)[0]
            except Exception as e:
                logger.warning("RF computation failed (mismatched leaf sets?): %s", e)
                rf = float('inf')
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
            node.add_face(TextFace(node.name, fsize=8), column=0, position="branch-right")

    def _build_a4_grid_layout(n_panels: int) -> tuple[list[dict], float, float]:
        # A4 portrait in points (72 dpi)
        page_w, page_h = 595.0, 842.0
        margin_x, margin_y = 20.0, 20.0
        gap_x, gap_y = 12.0, 12.0
        cols = 2
        rows = max(1, math.ceil(n_panels / cols))

        cell_w = (page_w - 2 * margin_x - gap_x) / cols
        cell_h = (page_h - 2 * margin_y - (rows - 1) * gap_y) / rows

        items = []
        for idx in range(n_panels):
            row = idx // cols
            col = idx % cols
            x0 = margin_x + col * (cell_w + gap_x)
            y0 = margin_y + row * (cell_h + gap_y)
            items.append(
                {"idx": idx, "x0": x0, "y0": y0, "x1": x0 + cell_w, "y1": y0 + cell_h}
            )
        return items, page_w, page_h

    def _merge_pdf_with_pymupdf(
        pdf_files: list[str],
        output_pdf_path: str,
        page_items: list[dict],
        page_width: float,
        page_height: float,
    ) -> bool:
        if fitz is None:
            return False
        merged_doc = fitz.open()
        page = merged_doc.new_page(width=page_width, height=page_height)
        for item in page_items:
            src = fitz.open(pdf_files[item["idx"]])
            src_rect = src[0].rect
            cell_w = item["x1"] - item["x0"]
            cell_h = item["y1"] - item["y0"]
            scale = min(cell_w / float(src_rect.width), cell_h / float(src_rect.height))
            draw_w = float(src_rect.width) * scale
            draw_h = float(src_rect.height) * scale
            dx = item["x0"] + (cell_w - draw_w) / 2.0
            dy = item["y0"] + (cell_h - draw_h) / 2.0
            rect = fitz.Rect(dx, dy, dx + draw_w, dy + draw_h)
            page.show_pdf_page(rect, src, 0)
            src.close()
        merged_doc.save(output_pdf_path)
        merged_doc.close()
        return True

    def _merge_pdf_with_pypdf(
        pdf_files: list[str],
        output_pdf_path: str,
        page_items: list[dict],
        page_width: float,
        page_height: float,
    ) -> bool:
        if PdfReader is None or PdfWriter is None or Transformation is None:
            return False
        writer = PdfWriter()
        page = writer.add_blank_page(width=page_width, height=page_height)
        for item in page_items:
            with open(pdf_files[item["idx"]], "rb") as handle:
                src_page = PdfReader(handle).pages[0]
                src_w = float(src_page.mediabox.width)
                src_h = float(src_page.mediabox.height)
                cell_w = item["x1"] - item["x0"]
                cell_h = item["y1"] - item["y0"]
                scale = min(cell_w / src_w, cell_h / src_h)
                draw_w = src_w * scale
                draw_h = src_h * scale
                tx = item["x0"] + (cell_w - draw_w) / 2.0
                # page_items y0/y1 are defined in a top-origin layout.
                # pypdf uses bottom-origin coordinates, so flip Y here to keep
                # panel order consistent with PyMuPDF (top-left -> bottom-right).
                ty = page_height - item["y0"] - draw_h - (cell_h - draw_h) / 2.0
                transform = Transformation().scale(scale, scale).translate(tx=tx, ty=ty)
                page.merge_transformed_page(src_page, transform)
        with open(output_pdf_path, "wb") as out_handle:
            writer.write(out_handle)
        return True

    if TreeStyle is None or NodeStyle is None or TextFace is None:
        logger.warning(
            "ete3 visualization components (TreeStyle/NodeStyle/TextFace) not available; "
            "skipping tree visualization."
        )
        return

    # Enforce descending frequency order for visualization panels.
    sorted_trees = sorted(
        tree_count_dict.items(),
        key=lambda item: (len(item[1]), item[0]),
        reverse=True,
    )[:top_n]
    temp_dir = tempfile.mkdtemp(prefix="temp_trees_")
    tree_panels_pdf = []
    for i, (tree_str, tree_group) in enumerate(sorted_trees):
        tree0 = read_tree(tree_str)
        tree = rename_input_tre(tree0, voucher2taxa_dic)
        tree.ladderize()
        tree.resolve_polytomy(recursive=True)
        tree.sort_descendants("support")

        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0
        nstyle["hz_line_type"] = 0
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        for node in tree.traverse():
            node.set_style(nstyle)

        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.layout_fn = _layout_fixed_leaf_font
        ts.title.add_face(TextFace(f"Topology {i} Gene trees: {len(tree_group)}", fsize=8), column=1)
        ts.show_scale = False
        # Keep branch-length geometry so ultrametric tip alignment is visible.
        ts.force_topology = False
        ts.draw_guiding_lines = True
        ts.guiding_lines_type = 2
        ts.guiding_lines_color = "#666666"
        ts.branch_vertical_margin = 8
        ts.margin_left = 24
        ts.margin_right = 24
        ts.margin_top = 20
        ts.margin_bottom = 20

        panel_pdf = os.path.join(temp_dir, f"tree_{i + 1}.pdf")
        tree.convert_to_ultrametric(tree_length=1)
        tree.render(panel_pdf, w=900, tree_style=ts)
        tree_panels_pdf.append(panel_pdf)

    if tree_panels_pdf:
        output_pdf = output_path if output_path.lower().endswith(".pdf") else f"{os.path.splitext(output_path)[0]}.pdf"
        page_items, page_w, page_h = _build_a4_grid_layout(len(tree_panels_pdf))
        merged_ok = _merge_pdf_with_pymupdf(tree_panels_pdf, output_pdf, page_items, page_w, page_h)
        if not merged_ok:
            merged_ok = _merge_pdf_with_pypdf(tree_panels_pdf, output_pdf, page_items, page_w, page_h)
        if not merged_ok:
            raise RuntimeError("Failed to merge vector PDF panels: install pymupdf or pypdf.")

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
            logger.warning(
                "Skipping tree %s: not single-copy (found %d copies for %d species)",
                k, len(t2), len(get_species_set(t2))
            )
            continue
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
