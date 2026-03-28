"""
Gene duplication loss visualization for the PhyloTracer pipeline.

This module summarizes duplication losses from tabular outputs and renders
annotated species trees with loss statistics and legends.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict

logger = logging.getLogger(__name__)

try:
    from ete3 import NodeStyle, TextFace, TreeStyle
except ImportError:
    NodeStyle = None
    TextFace = None
    TreeStyle = None


def identify_loss_detail(path_str: str):
    """Parse a loss path string and identify path-level loss transitions."""
    if not path_str or path_str == "NA":
        return []

    parsed_steps = []
    for step in path_str.split("->"):
        match = re.search(r"(.+?)\((\d+)\)$", step.strip())
        if match:
            parsed_steps.append((match.group(1).strip(), int(match.group(2))))

    if len(parsed_steps) < 2:
        return []

    losses = []
    for idx in range(1, len(parsed_steps)):
        prev_node, prev_copy = parsed_steps[idx - 1]
        curr_node, curr_copy = parsed_steps[idx]
        if curr_copy >= prev_copy:
            continue
        if prev_copy == 2 and curr_copy == 0:
            loss_type = "2-0"
        elif prev_copy == 2 and curr_copy == 1:
            loss_type = "2-1"
        elif prev_copy == 1 and curr_copy == 0:
            loss_type = "1-0"
        else:
            loss_type = "2-0" if (prev_copy - curr_copy) >= 2 else "2-1"
        losses.append((curr_node, loss_type))
    return losses


def parse_node_loss_events(node_events_str: str):
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


def _read_summary_rows(filepath: str):
    with open(filepath, 'r') as handle:
        header = handle.readline().strip().split('	')
        header_idx = {name: idx for idx, name in enumerate(header)}
        rows = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            cols = line.split('	')
            rows.append((cols, header_idx))
    return rows


def infer_visualizer_mode(filepath: str) -> str:
    """Infer how node-loss counts should be described in the figure.

    GD_Loss_Visualizer renders GD_Loss_Tracker output directly and does not
    ask the user to restate the counting mode. In practice, the tracker
    default is parsimony, so files without an explicit accumulate hint are
    interpreted as parsimony-style node counts.
    """
    if 'accumulate' in str(filepath).lower():
        return 'accumulate'
    return 'parsimony'


def get_event_stats(filepath: str):
    """Compute GD births and event-level 2-2/2-1/2-0 counts by GD birth node."""
    gd_birth_sets = defaultdict(set)
    best_event_type_by_level = defaultdict(dict)
    priority = {"2-2": 0, "2-1": 1, "2-0": 2}

    logger.info('Reading file: %s...', filepath)
    for cols, idx in _read_summary_rows(filepath):
        if len(cols) < 4:
            continue
        tree_id = cols[0]
        gd_id = cols[1]
        level_node = cols[3]
        event_key = f"{tree_id}|{gd_id}"
        gd_birth_sets[level_node].add(event_key)
        loss_type_i = idx.get('loss_type')
        if loss_type_i is None or loss_type_i >= len(cols):
            continue
        event_type = cols[loss_type_i].strip()
        if event_type not in priority:
            continue
        prev = best_event_type_by_level[level_node].get(event_key)
        if prev is None or priority[event_type] > priority[prev]:
            best_event_type_by_level[level_node][event_key] = event_type

    gd_births = {node: len(keys) for node, keys in gd_birth_sets.items()}
    event_loss_data = {}
    for node, event_map in best_event_type_by_level.items():
        counts = {'2-0': 0, '2-1': 0, '2-2': 0}
        for event_type in event_map.values():
            counts[event_type] += 1
        event_loss_data[node] = counts
    return gd_births, event_loss_data


def get_cumulative_counts_from_paths(filepath: str):
    """Aggregate cumulative path-level losses from ``loss_path``."""
    counts = defaultdict(lambda: {'2-0': 0, '2-1': 0, '1-0': 0, 'total': 0})
    for cols, idx in _read_summary_rows(filepath):
        path_i = idx.get('loss_path', 7)
        if path_i >= len(cols):
            continue
        for node_name, loss_type in identify_loss_detail(cols[path_i].strip()):
            counts[node_name][loss_type] += 1
            counts[node_name]['total'] += 1
    return dict(counts)


def get_node_event_counts(filepath: str):
    """Aggregate node-level path_count_node_events directly from tracker output."""
    counts = defaultdict(lambda: {'2-0': 0, '2-1': 0, '1-0': 0, 'total': 0})
    used = False
    for cols, idx in _read_summary_rows(filepath):
        events_i = idx.get('path_count_node_events')
        if events_i is None or events_i >= len(cols):
            continue
        events = parse_node_loss_events(cols[events_i].strip())
        if events:
            used = True
        for node_name, loss_type in events:
            counts[node_name][loss_type] += 1
            counts[node_name]['total'] += 1
    return dict(counts), used


def _get_maximal_cover_nodes(start_node, target_species):
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


def get_fallback_parsimony_counts(filepath: str, sptree):
    """Fallback parsimony-like compression from final loss_type assignments."""
    event_rows = {}
    for cols, idx in _read_summary_rows(filepath):
        loss_type_i = idx.get('loss_type')
        species_i = idx.get('species', 4)
        level_i = idx.get('level', 3)
        if loss_type_i is None or max(loss_type_i, species_i, level_i) >= len(cols):
            continue
        event_key = (cols[0], cols[1])
        slot = event_rows.setdefault(event_key, {
            'level': cols[level_i],
            'loss_one_species': set(),
            'loss_two_species': set(),
        })
        loss_type = cols[loss_type_i].strip()
        species = cols[species_i]
        if loss_type == '2-0':
            slot['loss_two_species'].add(species)
        elif loss_type == '2-1':
            slot['loss_one_species'].add(species)

    parsimony_counts = defaultdict(lambda: {'2-0': 0, '2-1': 0, '1-0': 0, 'total': 0})
    for event in event_rows.values():
        try:
            gd_node = sptree & event['level']
        except Exception:
            continue
        loss_two_cover = _get_maximal_cover_nodes(gd_node, event['loss_two_species'])
        covered_by_two = set()
        for node_name in loss_two_cover:
            parsimony_counts[node_name]['1-0'] += 1
            parsimony_counts[node_name]['total'] += 1
            try:
                covered_by_two.update((sptree & node_name).get_leaf_names())
            except Exception:
                pass
        residual_loss_one = set(event['loss_one_species']) - covered_by_two
        for node_name in _get_maximal_cover_nodes(gd_node, residual_loss_one):
            parsimony_counts[node_name]['2-1'] += 1
            parsimony_counts[node_name]['total'] += 1
    return dict(parsimony_counts)


def calculate_node_loss_score(event_loss_stats: dict, gd_birth_count: int):
    severe = event_loss_stats.get('2-0', 0)
    partial = event_loss_stats.get('2-1', 0)
    raw_score = (1.0 * severe) + (0.5 * partial)
    if gd_birth_count > 0:
        return raw_score, raw_score / gd_birth_count
    return raw_score, None


def format_node_loss_breakdown(node_stats: dict, loss_label: str) -> str:
    """Format per-node loss-type breakdown for display on the tree."""
    return (
        f"{loss_label}20/21/10="
        f"{node_stats.get('2-0', 0)}/"
        f"{node_stats.get('2-1', 0)}/"
        f"{node_stats.get('1-0', 0)}"
    )


def visualizer_sptree(
    filepath,
    sptree,
    output_file='gd_loss_pie_visualizer.PDF',
):
    """Render GD loss statistics on a species tree from one tracker summary.

    The visualizer consumes one GD_Loss_Tracker result at a time. Node-loss counts
    are taken directly from that summary. By default the summary is interpreted
    with parsimony semantics; accumulate-style wording is only used when the
    input path explicitly indicates accumulate output.
    """
    try:
        from ete3 import NodeStyle, TextFace, TreeStyle
    except ImportError:
        logger.error('ete3 is required for visualization. Install with: pip install ete3')
        return

    gd_births, event_loss_data = get_event_stats(filepath)

    visualizer_mode = infer_visualizer_mode(filepath)
    node_loss_counts, has_node_events = get_node_event_counts(filepath)
    if not has_node_events:
        if visualizer_mode == 'accumulate':
            node_loss_counts = get_cumulative_counts_from_paths(filepath)
        else:
            node_loss_counts = get_fallback_parsimony_counts(filepath, sptree)

    if visualizer_mode == 'accumulate':
        loss_label = 'C'
        loss_desc = 'accumulated descendant loss-path count from GD_Loss_Tracker'
        face_color = 'darkorange'
    else:
        loss_label = 'P'
        loss_desc = 'parsimony node-loss count from GD_Loss_Tracker (shared MRCA losses counted once)'
        face_color = 'darkgreen'

    logger.info('Node-loss source: %s', filepath)
    logger.info('Visualizer mode: %s', visualizer_mode)

    sptree.ladderize()
    sptree.sort_descendants('support')

    logger.info(
        '\n%-15s | %-8s | %-8s | %-8s | %-10s | %-10s | %-8s | %-8s',
        'Node', 'Birth', 'E2-0', 'E2-1', 'NodeLoss', 'Label', 'Lraw', 'Lnorm'
    )
    logger.info('-' * 92)

    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle['size'] = 0
        node.set_style(nstyle)
        node_name = node.name

        if node_name in gd_births and gd_births[node_name] > 0:
            node.add_face(
                TextFace(f"B={gd_births[node_name]}", fsize=6, fgcolor='blue'),
                column=0,
                position='branch-top',
            )

        if node_name in gd_births or node_name in event_loss_data or node_name in node_loss_counts:
            event_stats = event_loss_data.get(node_name, {'2-0': 0, '2-1': 0, '2-2': 0})
            node_stats = node_loss_counts.get(node_name, {'2-0': 0, '2-1': 0, '1-0': 0, 'total': 0})
            gd_birth_count = gd_births.get(node_name, 0)
            raw_loss_score, norm_loss_score = calculate_node_loss_score(event_stats, gd_birth_count)
            e_2_0 = event_stats.get('2-0', 0)
            e_2_1 = event_stats.get('2-1', 0)
            e_2_2 = event_stats.get('2-2', 0)
            node_total = node_stats.get('total', 0)

            node.add_face(TextFace(f"E={e_2_2}/{e_2_1}/{e_2_0}", fsize=6, fgcolor='red'), column=0, position='branch-bottom')
            if norm_loss_score is not None:
                node.add_face(TextFace(f"L={norm_loss_score:.2f}", fsize=6, fgcolor='purple'), column=0, position='branch-bottom')
            else:
                node.add_face(TextFace(f"Lraw={raw_loss_score:.2f}", fsize=6, fgcolor='gray'), column=0, position='branch-bottom')
            if node_total > 0:
                node.add_face(TextFace(f"{loss_label}={node_total}", fsize=6, fgcolor=face_color), column=0, position='branch-bottom')
                node.add_face(
                    TextFace(format_node_loss_breakdown(node_stats, loss_label), fsize=6, fgcolor=face_color),
                    column=0,
                    position='branch-bottom',
                )

            if node_total > 0 or e_2_0 > 0 or e_2_1 > 0:
                logger.info(
                    '%-15s | %-8d | %-8d | %-8d | %-10d | %-10s | %-8.2f | %-8.2f',
                    node_name, gd_birth_count, e_2_0, e_2_1, node_total, loss_label,
                    raw_loss_score, norm_loss_score if norm_loss_score is not None else float('nan')
                )

    ts = TreeStyle()
    ts.scale = 40
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.show_border = False
    ts.extra_branch_line_type = 0
    ts.legend_position = 1
    ts.title.add_face(TextFace('Legend:', bold=True, fsize=10), column=0)
    ts.title.add_face(TextFace('  B=Birth: GD events originating at this node', fsize=6, fgcolor='blue'), column=0)
    ts.title.add_face(TextFace('  Red E=2-2/2-1/2-0: final event-level copy-state counts', fsize=6, fgcolor='red'), column=0)
    ts.title.add_face(TextFace('  Purple L: event-level loss score = (E2-0 + 0.5*E2-1) / Birth', fsize=6, fgcolor='purple'), column=0)
    ts.title.add_face(TextFace('  Gray Lraw: shown when Birth=0 on a node', fsize=6, fgcolor='gray'), column=0)
    ts.title.add_face(TextFace(f'  {loss_label}: {loss_desc}', fsize=6, fgcolor=face_color), column=0)
    ts.title.add_face(TextFace(f'  {loss_label}20/21/10: node-loss breakdown by 2-0, 2-1, and 1-0 counts', fsize=6, fgcolor=face_color), column=0)
    ts.title.add_face(TextFace(f'  Note: Birth and {loss_label} quantify different processes and are not expected to match', fsize=6, fgcolor='black'), column=0)

    try:
        sptree.convert_to_ultrametric()
    except Exception:
        pass

    for leaf in sptree.iter_leaves():
        leaf.add_face(TextFace(leaf.name, fsize=7, fgcolor='black', fstyle='italic'), column=0)

    sptree.render(output_file, w=260, units='mm', tree_style=ts)
    logger.info('Visualization saved to %s', output_file)


if __name__ == '__main__':
    import argparse

    from ete3 import Tree

    parser = argparse.ArgumentParser(description='Visualize gene duplication losses on a species tree.')
    parser.add_argument('loss_summary_file', help='Path to one GD loss summary file.')
    parser.add_argument('species_tree_file', help='Path to the species tree file (Newick format).')
    parser.add_argument('-o', '--output', default='gd_loss_pie_visualizer.PDF', help='Output PDF file path (default: gd_loss_pie_visualizer.PDF).')
    args = parser.parse_args()

    sptree = Tree(args.species_tree_file, format=1)
    visualizer_sptree(
        args.loss_summary_file,
        sptree,
        output_file=args.output,
    )
