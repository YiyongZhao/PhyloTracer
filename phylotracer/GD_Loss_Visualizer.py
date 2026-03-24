"""
Gene duplication loss visualization for the PhyloTracer pipeline.

This module summarizes duplication losses from tabular outputs and renders
annotated species trees with loss statistics and legends.
"""

import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

try:
    from ete3 import CircleFace, NodeStyle, PieChartFace, TextFace, TreeStyle
except ImportError:
    CircleFace = None
    NodeStyle = None
    PieChartFace = None
    TextFace = None
    TreeStyle = None


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
        List of ``(node, loss_type)`` tuples.
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
            if prev_copy == 2 and curr_copy == 0:
                loss_type = "2-0"
            elif prev_copy == 2 and curr_copy == 1:
                loss_type = "2-1"
            elif prev_copy == 1 and curr_copy == 0:
                loss_type = "1-0"
            else:
                diff = prev_copy - curr_copy
                loss_type = "2-0" if diff >= 2 else "2-1"
            identified_losses.append((curr_node, loss_type))

    return identified_losses


def parse_node_loss_events(node_events_str):
    """Parse node-aware loss events encoded as ``node|type;node|type``."""
    if not node_events_str or node_events_str == "NA":
        return []
    events = []
    for item in node_events_str.split(";"):
        item = item.strip()
        if not item or "|" not in item:
            continue
        node_name, loss_type = item.split("|", 1)
        node_name = node_name.strip()
        loss_type = loss_type.strip()
        if node_name and loss_type in {"2-0", "2-1", "1-0"}:
            events.append((node_name, loss_type))
    return events


def get_stats_deduplicated(filepath):
    """
    Deduplicate GD events and compute birth and cumulative loss statistics.

    Returns
    -------
    tuple
        ``(gd_births, cumulative_losses, event_loss_counts)``.
    """
    gd_birth_sets = defaultdict(set)
    cumulative_loss_tracker = defaultdict(lambda: {"2-0": 0, "2-1": 0, "1-0": 0})
    event_loss_tracker = defaultdict(lambda: {"2-0": set(), "2-1": set(), "2-2": set()})
    event_type_priority = {"2-2": 0, "2-1": 1, "2-0": 2}
    best_event_type_by_level = defaultdict(dict)

    logger.info("Reading file: %s...", filepath)
    with open(filepath, "r") as f:
        header = f.readline().strip()
        if "tree_ID" not in header:
            raise ValueError(
                "GD_Loss_Visualizer currently requires gd_loss_summary.txt "
                "(detailed table with tree_ID column)."
            )
        header_cols = header.split("\t")
        header_idx = {name: idx for idx, name in enumerate(header_cols)}
        loss_type_idx = header_idx.get("loss_type")
        loss_path_idx = header_idx.get("loss_path", 7)

        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) <= loss_path_idx:
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
                    if prev is None or event_type_priority[event_loss_type] > event_type_priority[prev]:
                        best_event_type_by_level[level_node][unique_key] = event_loss_type

            for loss_node, loss_type in identify_loss_detail(cols[loss_path_idx].strip()):
                cumulative_loss_tracker[loss_node][loss_type] += 1

    final_births = {k: len(v) for k, v in gd_birth_sets.items()}
    final_cumulative_losses = {
        node: {
            "2-0": stats["2-0"],
            "2-1": stats["2-1"],
            "1-0": stats["1-0"],
        }
        for node, stats in cumulative_loss_tracker.items()
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

    return final_births, final_cumulative_losses, final_event_losses


def _get_maximal_cover_nodes(start_node, target_species):
    """Return maximal clades under ``start_node`` fully covered by ``target_species``."""
    target_species = set(target_species)
    if not target_species:
        return []

    result = []

    def _walk(node):
        leaves = set(node.get_leaf_names())
        if not (leaves & target_species):
            return
        if leaves.issubset(target_species):
            result.append(node.name)
            return
        for child in node.children:
            _walk(child)

    _walk(start_node)
    return result


def get_parsimony_loss_counts(filepath, sptree):
    """
    Compute simplified parsimony loss placements on the species tree.

    The rule is intentionally conservative for visualization:
    - species with final ``2-0`` are compressed to maximal shared clades and count as one
      parsimony loss placement per clade;
    - species with final ``2-1`` are compressed similarly, after removing leaves already
      explained by a ``2-0`` parsimony placement.

    This avoids repeated counting for situations such as A/B/C all losing in the same
    ancestral clade, where the figure should show a single ancestral loss placement.
    """
    event_rows = {}
    with open(filepath, "r") as f:
        header = f.readline().strip().split("\t")
        header_idx = {name: idx for idx, name in enumerate(header)}
        loss_type_idx = header_idx.get("loss_type")
        species_idx = header_idx.get("species", 4)
        level_idx = header_idx.get("level", 3)

        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if loss_type_idx is None or len(cols) <= loss_type_idx:
                continue
            event_key = (cols[0], cols[1])
            level = cols[level_idx]
            species = cols[species_idx]
            loss_type = cols[loss_type_idx].strip()
            slot = event_rows.setdefault(
                event_key,
                {"level": level, "loss_one_species": set(), "loss_two_species": set()},
            )
            if loss_type == "2-0":
                slot["loss_two_species"].add(species)
            elif loss_type == "2-1":
                slot["loss_one_species"].add(species)

    parsimony_counts = defaultdict(lambda: {"loss_one": 0, "loss_two": 0, "total": 0})
    for _, event in event_rows.items():
        level = event["level"]
        try:
            gd_node = sptree & level
        except Exception:
            continue

        loss_two_cover = _get_maximal_cover_nodes(gd_node, event["loss_two_species"])
        covered_by_two = set()
        for node_name in loss_two_cover:
            parsimony_counts[node_name]["loss_two"] += 1
            parsimony_counts[node_name]["total"] += 1
            try:
                covered_by_two.update((sptree & node_name).get_leaf_names())
            except Exception:
                pass

        residual_loss_one = set(event["loss_one_species"]) - covered_by_two
        loss_one_cover = _get_maximal_cover_nodes(gd_node, residual_loss_one)
        for node_name in loss_one_cover:
            parsimony_counts[node_name]["loss_one"] += 1
            parsimony_counts[node_name]["total"] += 1

    return dict(parsimony_counts)


def calculate_node_loss_score(event_loss_stats: dict, gd_birth_count: int):
    """Compute node-level raw and normalized GD-loss severity scores."""
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
    """Print loss statistics along the path to a target species."""
    target_nodes = sptree.search_nodes(name=target_species)
    if not target_nodes:
        logger.error("Species '%s' not found in tree", target_species)
        return

    leaf = target_nodes[0]
    path_nodes = leaf.get_ancestors()
    path_nodes.reverse()
    path_nodes.append(leaf)

    logger.info("\n%s Path: Root -> %s %s", '=' * 20, target_species, '=' * 20)
    logger.info("%-20s | %-10s | %-10s | %-10s | %s", 'Node', '2->0', '2->1', '1->0', 'Total')
    logger.info("-" * 80)

    for node in path_nodes:
        node_name = node.name
        if node_name in final_losses:
            stats = final_losses[node_name]
            m = stats["2-0"]
            n1 = stats["2-1"]
            n0 = stats["1-0"]
            total = m + n1 + n0
            if total > 0:
                logger.info("%-20s | %-10d | %-10d | %-10d | %d", node_name, m, n1, n0, total)
    logger.info("=" * 75)


# ======================================================
# Section 3: Visualization Rendering
# ======================================================


def visualizer_sptree(filepath, sptree, output_file="gd_loss_pie_visualizer.PDF"):
    """Render GD loss statistics on a species tree."""
    try:
        from ete3 import TreeStyle, NodeStyle
    except ImportError:
        logger.error("ete3 is required for visualization. Install with: pip install ete3")
        return

    gd_births, cumulative_loss_data, event_loss_data = get_stats_deduplicated(filepath)
    parsimony_loss_data = get_parsimony_loss_counts(filepath, sptree)

    sptree.ladderize()
    sptree.sort_descendants("support")

    logger.info(
        "\n%-15s | %-8s | %-8s | %-8s | %-10s | %-10s | %-8s | %-8s",
        "Node",
        "Birth",
        "E2-0",
        "E2-1",
        "CumLoss",
        "ParsLoss",
        "Lraw",
        "Lnorm",
    )
    logger.info("-" * 92)

    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        node.set_style(nstyle)
        node_name = node.name

        if node_name in gd_births and gd_births[node_name] > 0:
            node.add_face(TextFace(f"{gd_births[node_name]}", fsize=6, fgcolor="blue"), column=0, position="branch-top")

        if node_name in cumulative_loss_data or node_name in gd_births or node_name in event_loss_data or node_name in parsimony_loss_data:
            cumulative_stats = cumulative_loss_data.get(node_name, {"2-0": 0, "2-1": 0, "1-0": 0})
            c_total = cumulative_stats["2-0"] + cumulative_stats["2-1"] + cumulative_stats["1-0"]
            parsimony_stats = parsimony_loss_data.get(node_name, {"loss_one": 0, "loss_two": 0, "total": 0})
            p_total = parsimony_stats["total"]
            gd_birth_count = gd_births.get(node_name, 0)
            event_stats = event_loss_data.get(node_name, {"2-0": 0, "2-1": 0, "2-2": 0})
            raw_loss_score, norm_loss_score = calculate_node_loss_score(event_stats, gd_birth_count)
            e_2_0 = event_stats.get("2-0", 0)
            e_2_1 = event_stats.get("2-1", 0)
            e_2_2 = event_stats.get("2-2", 0)

            node.add_face(TextFace(f"E={e_2_2}/{e_2_1}/{e_2_0}", fsize=6, fgcolor="red"), column=0, position="branch-bottom")
            if norm_loss_score is not None:
                node.add_face(TextFace(f"L={norm_loss_score:.2f}", fsize=6, fgcolor="purple"), column=0, position="branch-bottom")
            else:
                node.add_face(TextFace(f"Lraw={raw_loss_score:.2f}", fsize=6, fgcolor="gray"), column=0, position="branch-bottom")
            if c_total > 0:
                node.add_face(TextFace(f"C={c_total}", fsize=6, fgcolor="darkorange"), column=0, position="branch-bottom")
            if p_total > 0:
                node.add_face(TextFace(f"P={p_total}", fsize=6, fgcolor="darkgreen"), column=0, position="branch-bottom")

            if c_total > 0 or p_total > 0 or e_2_0 > 0 or e_2_1 > 0:
                logger.info(
                    "%-15s | %-8d | %-8d | %-8d | %-10d | %-10d | %-8.2f | %-8.2f",
                    node_name,
                    gd_birth_count,
                    e_2_0,
                    e_2_1,
                    c_total,
                    p_total,
                    raw_loss_score,
                    norm_loss_score if norm_loss_score is not None else float("nan"),
                )

    ts = TreeStyle()
    ts.scale = 40
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.show_border = False
    ts.extra_branch_line_type = 0
    ts.legend_position = 1
    ts.title.add_face(TextFace("Legend:", bold=True, fsize=10), column=0)
    ts.title.add_face(TextFace("  Blue number (top): GD events generated on this node", fsize=6, fgcolor="blue"), column=0)
    ts.title.add_face(TextFace("  Red E=2-2/2-1/2-0: event-level final copy-state counts", fsize=6, fgcolor="red"), column=0)
    ts.title.add_face(TextFace("  Purple L: event-level loss score = (E2-0 + 0.5*E2-1) / GD_birth", fsize=6, fgcolor="purple"), column=0)
    ts.title.add_face(TextFace("  Gray Lraw: shown when GD_birth=0 on a node", fsize=6, fgcolor="gray"), column=0)
    ts.title.add_face(TextFace("  Orange C: cumulative path-loss count (all descendant paths counted)", fsize=6, fgcolor="darkorange"), column=0)
    ts.title.add_face(TextFace("  Green P: parsimony loss count (shared descendant losses compressed to ancestor clades)", fsize=6, fgcolor="darkgreen"), column=0)

    try:
        sptree.convert_to_ultrametric()
    except Exception as exc:
        logger.debug("convert_to_ultrametric failed (non-fatal): %s", exc)

    for leaf in sptree.iter_leaves():
        leaf.add_face(TextFace(leaf.name, fsize=7, fgcolor="black", fstyle="italic"), column=0)

    sptree.render(output_file, w=260, units="mm", tree_style=ts)
    logger.info("Visualization saved to %s", output_file)


# ======================================================
# Section 4: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse
    from ete3 import Tree

    parser = argparse.ArgumentParser(description="Visualize gene duplication losses on a species tree.")
    parser.add_argument("loss_summary_file", help="Path to the loss summary file (tabular format).")
    parser.add_argument("species_tree_file", help="Path to the species tree file (Newick format).")
    parser.add_argument("-o", "--output", default="gd_loss_pie_visualizer.PDF", help="Output PDF file path (default: gd_loss_pie_visualizer.PDF).")
    args = parser.parse_args()

    sptree = Tree(args.species_tree_file)
    visualizer_sptree(args.loss_summary_file, sptree, output_file=args.output)
