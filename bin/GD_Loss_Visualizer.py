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
    event_loss_tracker = defaultdict(lambda: {"2-0": set(), "2-1": set(), "2-2": set()})
    # Enforce mutually exclusive event-level type by keeping the most severe type
    # per (level, event_key): 2-0 > 2-1 > 2-2.
    event_type_priority = {"2-2": 0, "2-1": 1, "2-0": 2}
    best_event_type_by_level = defaultdict(dict)
    print(f"Reading file: {filepath}...")
    with open(filepath, "r") as f:
        header = f.readline().strip()
        if "tree_ID" not in header:
            raise ValueError(
                "GD_Loss_Visualizer currently requires gd_loss_summary.txt "
                "(detailed table with tree_ID column)."
            )
        header_cols = header.split("\t")
        header_idx = {name: idx for idx, name in enumerate(header_cols)}
        loss_type_idx = header_idx.get("loss_type", None)

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
            if loss_type_idx is not None and loss_type_idx < len(cols):
                event_loss_type = cols[loss_type_idx].strip()
                if event_loss_type in event_type_priority:
                    prev = best_event_type_by_level[level_node].get(unique_key)
                    if prev is None:
                        best_event_type_by_level[level_node][unique_key] = event_loss_type
                    elif event_type_priority[event_loss_type] > event_type_priority[prev]:
                        best_event_type_by_level[level_node][unique_key] = event_loss_type

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

    final_event_losses = {}
    for node, event_type_map in best_event_type_by_level.items():
        for event_key, event_type in event_type_map.items():
            event_loss_tracker[node][event_type].add(event_key)
        final_event_losses[node] = {
            "2-0": len(event_loss_tracker[node]["2-0"]),
            "2-1": len(event_loss_tracker[node]["2-1"]),
            "2-2": len(event_loss_tracker[node]["2-2"]),
        }

    return final_births, final_losses, final_event_losses


def calculate_node_loss_score(event_loss_stats: dict, gd_birth_count: int):
    """Compute node-level raw and normalized GD-loss severity scores.

    The raw score emphasizes severe loss more than partial loss:
    ``raw = 1*(2->0) + 0.5*(2->1)``.
    The normalized score is ``raw / GD_birth`` when ``GD_birth > 0``.
    """
    severe = event_loss_stats.get("2-0", 0)
    partial = event_loss_stats.get("2-1", 0)
    raw_score = (1.0 * severe) + (0.5 * partial)
    if gd_birth_count > 0:
        return raw_score, raw_score / gd_birth_count
    return raw_score, None


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


def visualizer_sptree(
    filepath,
    sptree,
    output_file="gd_loss_pie_visualizer.PDF",
):
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
    gd_births, loss_data, event_loss_data = get_stats_deduplicated(filepath)

    sptree.ladderize()
    sptree.sort_descendants("support")

    colors = {
        "2-0": "#D62728",
        "2-1": "#FF7F0E",
        "1-0": "#1F77B4",
    }
    color_list = [colors["2-0"], colors["2-1"], colors["1-0"]]

    print(
        f"\n{'Node':<15} | {'GD Birth':<10} | {'E2-0':<8} | {'E2-1':<8} | "
        f"{'Path2-0':<8} | {'Path2-1':<8} | {'Path1-0':<8} | "
        f"{'RawLoss':<8} | {'NormLoss':<8}"
    )
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

        if node_name in loss_data or node_name in gd_births or node_name in event_loss_data:
            stats = loss_data.get(node_name, {"2-0": 0, "2-1": 0, "1-0": 0})
            c_2_0 = stats["2-0"]
            c_2_1 = stats["2-1"]
            c_1_0 = stats["1-0"]
            gd_birth_count = gd_births.get(node_name, 0)
            event_stats = event_loss_data.get(node_name, {"2-0": 0, "2-1": 0, "2-2": 0})
            raw_loss_score, norm_loss_score = calculate_node_loss_score(
                event_stats, gd_birth_count
            )

            total_loss = c_2_0 + c_2_1 + c_1_0
            # Show event-level counts in red so they are consistent with L formula.
            e_2_0 = event_stats.get("2-0", 0)
            e_2_1 = event_stats.get("2-1", 0)
            e_2_2 = event_stats.get("2-2", 0)
            node.add_face(
                TextFace(f"{e_2_2}/{e_2_1}/{e_2_0}", fsize=6, fgcolor="red"),
                column=0,
                position="branch-bottom",
            )
            if norm_loss_score is not None:
                node.add_face(
                    TextFace(f"L={norm_loss_score:.2f}", fsize=6, fgcolor="purple"),
                    column=0,
                    position="branch-bottom",
                )
            else:
                node.add_face(
                    TextFace(f"Lraw={raw_loss_score:.2f}", fsize=6, fgcolor="gray"),
                    column=0,
                    position="branch-bottom",
                )
            if total_loss > 0:
                print(
                    f"{node_name:<15} | {gd_births.get(node_name,0):<10} | "
                    f"{event_stats.get('2-0', 0):<8} | {event_stats.get('2-1', 0):<8} | "
                    f"{c_2_0:<8} | {c_2_1:<8} | {c_1_0:<8} | "
                    f"{raw_loss_score:<8.2f} | "
                    f"{(norm_loss_score if norm_loss_score is not None else float('nan')):<8.2f}"
                )

            # percents = [
            #     (c_2_0 / total_loss) * 100,
            #     (c_2_1 / total_loss) * 100,
            #     (c_1_0 / total_loss) * 100,
            # ]
            # pie = PieChartFace(percents, width=6, height=6, colors=color_list)
            # node.add_face(pie, column=0, position="branch-bottom")

    ts = TreeStyle()
    ts.scale = 40
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.show_border = False
    ts.extra_branch_line_type = 0
    

    ts.legend_position = 1
    ts.title.add_face(TextFace("Legend:", bold=True, fsize=10), column=0)
    ts.title.add_face(
        TextFace("  Blue Number (Top): GD Events Generated", fsize=6, fgcolor="blue"),
        column=0,
    )
    ts.title.add_face(
        TextFace("  Red (Bottom): E(2-2)/E(2-1)/E(2-0) event-level counts", fsize=6),
        column=0,
    )
    ts.title.add_face(
        TextFace(
            "  Purple L score (event-level): (E(2-0) + 0.5*E(2-1)) / GD_birth",
            fsize=6,
            fgcolor="purple",
        ),
        column=0,
    )
    ts.title.add_face(
        TextFace("  Gray Lraw: shown when GD_birth=0 on a node", fsize=6),
        column=0,
    )
    try:
        realign_branch_length(sptree)
        rejust_root_dist(sptree)
    except Exception:
        pass

    for leaf in sptree.iter_leaves():
        leaf.add_face(
            TextFace(leaf.name, fsize=7, fgcolor="black", fstyle="italic"),
            column=0)

    sptree.render(output_file, w=260, units="mm", tree_style=ts)
    print(f"\nVisualization saved to {output_file}")


# ======================================================
# Section 4: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    pass
