"""
Gene duplication loss visualization for the PhyloTracer pipeline.

This module summarizes duplication losses from tabular outputs and renders
annotated species trees with loss statistics and legends.
"""

import re
from collections import defaultdict

from __init__ import *
from ete3 import CircleFace, PieChartFace

# ======================================================
# Section 1: Loss Parsing and Statistics
# ======================================================


def identify_loss_detail(path_str):
    """
    Parse a loss path string and identify loss nodes and types.

    Parameters
    ----------
    path_str : str
        Path string formatted as ``node(count)->node(count)``.

    Returns
    -------
    list
        List of (node, loss_type) tuples.

    Assumptions
    -----------
    Path strings encode copy numbers in parentheses.
    """
    if not path_str or path_str == "NA":
        return []

    steps = path_str.split("->")
    parsed_steps = []

    for step in steps:
        step = step.strip()
        match = re.search(r"(.+?)\((\d+)\)$", step)
        if not match:
            continue
        parsed_steps.append((match.group(1).strip(), int(match.group(2))))

    if len(parsed_steps) < 2:
        return []

    identified_losses = []

    for i in range(1, len(parsed_steps)):
        prev_node, prev_copy = parsed_steps[i - 1]
        curr_node, curr_copy = parsed_steps[i]

        if curr_copy < prev_copy:
            l_type = None
            if prev_copy == 2 and curr_copy == 0:
                l_type = "2-0"
            elif prev_copy == 2 and curr_copy == 1:
                l_type = "2-1"
            elif prev_copy == 1 and curr_copy == 0:
                l_type = "1-0"
            else:
                diff = prev_copy - curr_copy
                l_type = "2-0" if diff >= 2 else "2-1"

            if l_type:
                identified_losses.append((curr_node, l_type))

    return identified_losses


def get_stats_deduplicated(filepath):
    """
    Deduplicate GD events and compute birth and loss statistics.

    Parameters
    ----------
    filepath : str
        Path to the loss summary file.

    Returns
    -------
    tuple
        (gd_births, loss_counts) dictionaries.

    Assumptions
    -----------
    Tree ID and GD ID together form a unique event key.
    """
    gd_birth_sets = defaultdict(set)
    loss_tracker = defaultdict(lambda: {"2-0": set(), "2-1": set(), "1-0": set()})

    print(f"Reading file: {filepath}...")
    with open(filepath, "r") as f:
        header = f.readline()
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue

            tree_id = cols[0]
            gd_id_raw = cols[1]
            unique_key = f"{tree_id}|{gd_id_raw}"

            level_node = cols[3]
            gd_birth_sets[level_node].add(unique_key)

            loss_path = cols[7].strip()
            all_losses = identify_loss_detail(loss_path)

            for loss_node, loss_type in all_losses:
                loss_tracker[loss_node][loss_type].add(unique_key)

    final_births = {k: len(v) for k, v in gd_birth_sets.items()}
    final_losses = {}
    for node, types in loss_tracker.items():
        final_losses[node] = {
            "2-0": len(types["2-0"]),
            "2-1": len(types["2-1"]),
            "1-0": len(types["1-0"]),
        }

    return final_births, final_losses


# ======================================================
# Section 2: Reporting Utilities
# ======================================================


def print_path_stats(sptree, final_losses, target_species="Arabidopsis_thaliana"):
    """
    Print loss statistics along the path to a target species.

    Parameters
    ----------
    sptree : object
        Species tree object.
    final_losses : dict
        Loss counts by node and type.
    target_species : str, optional
        Target species name for path reporting.

    Returns
    -------
    None

    Assumptions
    -----------
    The species name exists in the tree.
    """
    target_nodes = sptree.search_nodes(name=target_species)
    if not target_nodes:
        print(f"Error: Species '{target_species}' not found in tree")
        return

    leaf = target_nodes[0]
    path_nodes = leaf.get_ancestors()
    path_nodes.reverse()
    path_nodes.append(leaf)

    print(f"\n{'='*20} Path: Root -> {target_species} {'='*20}")
    print(f"{'Node':<20} | {'2->0':<10} | {'2->1':<10} | {'1->0':<10} | {'Total'}")
    print("-" * 80)

    for node in path_nodes:
        node_name = node.name
        if node_name in final_losses:
            stats = final_losses[node_name]
            m = stats["2-0"]
            n1 = stats["2-1"]
            n0 = stats["1-0"]
            total = m + n1 + n0
            if total > 0:
                print(f"{node_name:<20} | {m:<10} | {n1:<10} | {n0:<10} | {total}")
    print("=" * 75)


# ======================================================
# Section 3: Visualization Rendering
# ======================================================


def visualizer_sptree(filepath, sptree, output_file="gd_loss_pie_visualizer.PDF"):
    """
    Render GD loss statistics on a species tree.

    Parameters
    ----------
    filepath : str
        Path to the loss summary file.
    sptree : object
        Species tree object to annotate.
    output_file : str, optional
        Output PDF file path.

    Returns
    -------
    None

    Assumptions
    -----------
    Loss summary file follows the expected tabular format.
    """
    gd_births, loss_data = get_stats_deduplicated(filepath)

    sptree.ladderize()

    colors = {
        "2-0": "#D62728",
        "2-1": "#FF7F0E",
        "1-0": "#1F77B4",
    }
    color_list = [colors["2-0"], colors["2-1"], colors["1-0"]]

    print(f"\n{'Node':<15} | {'GD Birth':<10} | {'2->0':<8} | {'2->1':<8} | {'1->0':<8}")
    print("-" * 80)

    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        node.set_style(nstyle)

        node_name = node.name

        if node_name in gd_births and gd_births[node_name] > 0:
            cnt = gd_births[node_name]
            face = TextFace(f"{cnt}", fsize=6, fgcolor="blue")
            node.add_face(face, column=0, position="branch-top")

        if node_name in loss_data:
            stats = loss_data[node_name]
            c_2_0 = stats["2-0"]
            c_2_1 = stats["2-1"]
            c_1_0 = stats["1-0"]

            total_loss = c_2_0 + c_2_1 + c_1_0
            node.add_face(
                TextFace(f"{c_2_0}/{c_2_1}/{c_1_0}", fsize=6, fgcolor="red"),
                column=0,
                position="branch-bottom",
            )
            if total_loss > 0:
                print(
                    f"{node_name:<15} | {gd_births.get(node_name,0):<10} | "
                    f"{c_2_0:<8} | {c_2_1:<8} | {c_1_0:<8}"
                )

            # percents = [
            #     (c_2_0 / total_loss) * 100,
            #     (c_2_1 / total_loss) * 100,
            #     (c_1_0 / total_loss) * 100,
            # ]
            # pie = PieChartFace(percents, width=6, height=6, colors=color_list)
            # node.add_face(pie, column=0, position="branch-bottom")

    ts = TreeStyle()
    ts.scale = 20
    ts.show_leaf_name = True
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.show_border = False

    ts.title.add_face(TextFace("Legend:", bold=True, fsize=10), column=0)
    ts.title.add_face(
        TextFace("  Blue Number (Top): GD Events Generated", fsize=6, fgcolor="blue"),
        column=0,
    )
    ts.title.add_face(
        TextFace("  Pie Chart (Bottom): Loss Types Distribution", fsize=6),
        column=0,
    )

    c1 = CircleFace(radius=4, color=colors["2-0"], style="circle")
    t1 = TextFace(" 2->0 (Severe Loss)", fsize=6)
    ts.title.add_face(c1, column=0)
    ts.title.add_face(t1, column=1)

    c2 = CircleFace(radius=4, color=colors["2-1"], style="circle")
    t2 = TextFace(" 2->1 (Partial Loss)", fsize=6)
    ts.title.add_face(c2, column=0)
    ts.title.add_face(t2, column=1)

    c3 = CircleFace(radius=4, color=colors["1-0"], style="circle")
    t3 = TextFace(" 1->0 (Final Loss)", fsize=6)
    ts.title.add_face(c3, column=0)
    ts.title.add_face(t3, column=1)

    try:
        realign_branch_length(sptree)
        rejust_root_dist(sptree)
    except Exception:
        pass

    sptree.render(output_file, w=210, units="mm", tree_style=ts)
    print(f"\nVisualization saved to {output_file}")


# ======================================================
# Section 4: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    pass
