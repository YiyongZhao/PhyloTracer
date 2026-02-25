"""
Gene duplication visualization utilities for the PhyloTracer pipeline.

This module parses duplication results, aggregates events by species-tree
nodes, and renders annotated species trees with duplication summaries.
"""

import os
import re
from collections import defaultdict

try:
    from ete3 import CircleFace, NodeStyle, PieChartFace, TextFace, TreeStyle
except ImportError:
    CircleFace = None
    NodeStyle = None
    PieChartFace = None
    TextFace = None
    TreeStyle = None

from phylotracer import (
    read_phylo_tree,
    read_and_return_dict,
)

# ======================================================
# Section 1: GD Result Parsing and Aggregation
# ======================================================


def process_gd_result(gd_file: str) -> list:
    """
    Parse gene duplication results produced by the pipeline.

    Parameters
    ----------
    gd_file : str
        Path to the gene duplication result file.

    Returns
    -------
    list
        List of (gd_id, level, gd_type) tuples.

    Assumptions
    -----------
    The input file is tab-delimited with a header line starting with '#'.
    """
    records = []
    with open(gd_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 10:
                continue

            gd_id = parts[1]
            level = parts[5]
            gd_type_raw = parts[9]
            # Backward compatibility for historical labels.
            legacy_to_current = {
                "ABAB": "AABB",
                "ABB": "AXBB",
                "AAB": "AABX",
            }
            gd_type = legacy_to_current.get(gd_type_raw, gd_type_raw)

            records.append((gd_id, level, gd_type))

    return records


def get_count_dic(gene_duplications: list) -> dict:
    """
    Count unique duplication events per level and GD type.

    Parameters
    ----------
    gene_duplications : list
        List of (gd_id, level, gd_type) tuples.

    Returns
    -------
    dict
        Nested dictionary mapping level to GD-type counts.

    Assumptions
    -----------
    Each GD identifier is unique per level.
    """
    level_type_count = defaultdict(lambda: defaultdict(int))
    seen_ids_per_level = defaultdict(set)

    for gd_id, level, gd_type in gene_duplications:
        if gd_id not in seen_ids_per_level[level]:
            level_type_count[level][gd_type] += 1
            seen_ids_per_level[level].add(gd_id)

    return dict(level_type_count)


# ======================================================
# Section 2: Species Tree Annotation and Rendering
# ======================================================


def mark_sptree(
    sptree: object,
    count_dic: dict,
    taxa: dict,
    output_pdf: str = "phylotracer_gd_visualizer.pdf",
) -> object:
    """
    Annotate a species tree with duplication summaries and taxa labels.

    Parameters
    ----------
    sptree : object
        Species tree object to annotate and render.
    count_dic : dict
        Duplication counts by species-tree node and GD type.
    taxa : dict
        Mapping from leaf names to taxa labels.

    Returns
    -------
    object
        Rendered output from ``Tree.render``.

    Assumptions
    -----------
    Species tree node names match those used in duplication summaries.
    """
    sptree.ladderize()
    sptree.sort_descendants("support")

    ts = TreeStyle()
    ts.scale = 10
    ts.legend_position = 1
    ts.show_leaf_name = False
    ts.guiding_lines_type = 0
    ts.guiding_lines_color = "black"
    ts.draw_guiding_lines = True
    ts.extra_branch_line_type = 0
    ts.extra_branch_line_color = "black"

    ts.legend.add_face(TextFace("  Legend:", fsize=6, ftype="Arial"), column=0)
    ts.legend.add_face(
        TextFace(
            "  Red numbers: Gene duplication events",
            fsize=6,
            fgcolor="red",
            ftype="Arial",
        ),
        column=0,
    )
    ts.legend.add_face(
        TextFace(
            "  Blue numbers: Node identifiers",
            fsize=6,
            fgcolor="blue",
            ftype="Arial",
        ),
        column=0,
    )

    type_colors = {
        "AABB": "#1f77b4",
        "AXBB": "#ff7f0e",
        "AABX": "#2ca02c",
        "Complex": "#d62728",
    }
    type_order = ["AABB", "AXBB", "AABX", "Complex"]

    ts.title.add_face(TextFace(" GD Events Distribution ", fsize=6, ftype="Arial"), column=0)
    ts.title.add_face(CircleFace(4, type_colors["AABB"]), column=1)
    ts.title.add_face(TextFace(" AABB", fsize=6, ftype="Arial"), column=2)
    ts.title.add_face(CircleFace(4, type_colors["AXBB"]), column=3)
    ts.title.add_face(TextFace(" AXBB", fsize=6), column=4)
    ts.title.add_face(CircleFace(4, type_colors["AABX"]), column=5)
    ts.title.add_face(TextFace(" AABX", fsize=6, ftype="Arial"), column=6)
    ts.title.add_face(CircleFace(4, type_colors["Complex"]), column=7)
    ts.title.add_face(TextFace(" Complex", fsize=6, ftype="Arial"), column=8)



    sptree = sptree.copy()
    for leaf in sptree:
        if leaf.name in taxa:
            leaf.name = taxa[leaf.name]

    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        node.set_style(nstyle)
        node_name = node.name

        pie_face = None
        if node_name in count_dic:
            counts = {t: count_dic[node_name].get(t, 0) for t in type_order}
            values = [counts[t] for t in type_order]
            colors = [type_colors[t] for t in type_order]
            total = sum(values)
            if total > 0:
                percentages = [round(v / total * 100, 1) for v in values]
                diff = 100.0 - sum(percentages)
                if percentages:
                    percentages[0] += diff
                pie_face = PieChartFace(percentages, width=6, height=6, colors=colors)

        gd_text_face = None
        if node_name in count_dic:
            total_gd = sum(count_dic[node_name].values())
            if total_gd > 0:
                gd_text_face = TextFace(f"{total_gd}", fsize=6, fgcolor="red", ftype="Arial")

        if node.is_leaf():
            if pie_face:
                node.add_face(pie_face, column=1, position="branch-top")
            if gd_text_face:
                node.add_face(gd_text_face, column=0, position="branch-top")
            name_face = TextFace(node_name, fsize=6, fgcolor="black", ftype="Arial", fstyle="italic")
            node.add_face(name_face, column=1, position="aligned")

        else:
            if gd_text_face:
                node.add_face(gd_text_face, column=0, position="branch-top")
            if pie_face:
                node.add_face(pie_face, column=1, position="branch-top")

            id_face = TextFace(node_name, fsize=6, fgcolor="blue", ftype="Arial")
            node.add_face(id_face, column=0, position="branch-bottom")

    sptree.convert_to_ultrametric()

    return sptree.render(output_pdf, w=210, units="mm", tree_style=ts)


# ======================================================
# Section 3: Main Pipeline (Orchestration)
# ======================================================


def gd_visualizer_main(sptree, gd_result, taxa):
    """
    Visualize gene duplication events on a species tree.

    Parameters
    ----------
    sptree : object
        Species tree object to be visualized.
    gd_result : str
        Path to the gene duplication result file.
    taxa : dict
        Mapping from leaf names to taxa labels.

    Returns
    -------
    None

    Assumptions
    -----------
    The GD result file uses canonical GD labels per event.
    """
    gds = process_gd_result(gd_result)
    count_dic = get_count_dic(gds)
    gd_base = os.path.splitext(os.path.basename(gd_result))[0]
    output_pdf = f"{gd_base}.pdf"
    mark_sptree(sptree, count_dic, taxa, output_pdf=output_pdf)


# ======================================================
# Section 4: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Visualize gene duplication events on a species tree.",
    )
    parser.add_argument(
        "species_tree_file",
        help="Path to the species tree file (Newick format).",
    )
    parser.add_argument(
        "gd_result_file",
        help="Path to the gene duplication result file.",
    )
    parser.add_argument(
        "taxa_file",
        help="Path to the taxa mapping file.",
    )
    args = parser.parse_args()

    sptree = read_phylo_tree(args.species_tree_file)
    taxa_data = read_and_return_dict(args.taxa_file)

    gd_visualizer_main(sptree, args.gd_result_file, taxa_data)
