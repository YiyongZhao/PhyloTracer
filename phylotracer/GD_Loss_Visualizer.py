"""
Gene duplication loss visualization for the PhyloTracer pipeline.

This module summarizes duplication losses from tabular outputs and renders
annotated species trees with loss statistics and legends.
"""

from __future__ import annotations

import logging
import os
import re
from collections import defaultdict
import csv
logger = logging.getLogger(__name__)

from phylotracer import parse_loss_path, is_internal_node

try:
    from ete3 import NodeStyle, TextFace, TreeStyle
except ImportError:
    NodeStyle = None
    TextFace = None
    TreeStyle = None


VALID_VISUALIZER_MODES = {"parsimony", "accumulate"}


def normalize_visualizer_mode(mode, default="parsimony"):
    """Normalize one visualizer mode name."""
    text = str(mode or "").strip().lower()
    if text in VALID_VISUALIZER_MODES:
        return text
    return default


def parse_path_count_summary(summary_str: str):
    """Parse one merged path-count summary field from tracker output."""
    if not summary_str or summary_str == "NA":
        return {"node_events": "NA", "types": "NA"}

    parts = [p.strip() for p in str(summary_str).split(";") if p.strip()]
    node_events = []
    types = []
    for part in parts:
        if ":" in part:
            node_name, loss_type = part.split(":", 1)
            node_name = node_name.strip()
            loss_type = loss_type.strip().replace("-", "_")
            if node_name and loss_type in {"2_0", "2_1", "1_0"}:
                node_events.append(f"{node_name}|{loss_type}")
                types.append(loss_type)
        else:
            loss_type = part.strip().replace("-", "_")
            if loss_type in {"2_0", "2_1", "1_0"}:
                types.append(loss_type)
    return {
        "node_events": ";".join(node_events) if node_events else "NA",
        "types": ";".join(types) if types else "NA",
    }


def identify_loss_detail(path_str: str):
    """Parse a loss path string and identify path-level loss transitions."""
    parsed_steps = parse_loss_path(path_str)

    if len(parsed_steps) < 2:
        return []

    losses = []
    for idx in range(1, len(parsed_steps)):
        prev_node, prev_copy = parsed_steps[idx - 1]
        curr_node, curr_copy = parsed_steps[idx]
        if curr_copy >= prev_copy:
            continue
        if prev_copy == 2 and curr_copy == 0:
            loss_type = "2_0"
        elif prev_copy == 2 and curr_copy == 1:
            loss_type = "2_1"
        elif prev_copy == 1 and curr_copy == 0:
            loss_type = "1_0"
        else:
            loss_type = "2_0" if (prev_copy - curr_copy) >= 2 else "2_1"
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
        loss_type = loss_type.strip().replace("-", "_")
        if node_name and loss_type in {"2_0", "2_1", "1_0"}:
            events.append((node_name, loss_type))
    return events


def _read_summary_rows(filepath: str):
    delimiter = "," if str(filepath).lower().endswith(".csv") else "\t"
    with open(filepath, 'r', newline='') as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader, [])
        header_idx = {name: idx for idx, name in enumerate(header)}
        rows = []
        for cols in reader:
            if not cols:
                continue
            rows.append((cols, header_idx))
    return rows


def _find_col(idx: dict, *names: str):
    for name in names:
        if name in idx:
            return idx[name]
    return None


def find_node_summary_file(loss_result_file: str):
    """Return the tracker-generated node summary path next to ``gd_loss.csv``."""
    candidate = os.path.join(
        os.path.dirname(os.path.abspath(loss_result_file)),
        "gd_loss_node_summary.tsv",
    )
    return candidate if os.path.exists(candidate) else None


def load_node_summary_counts(node_summary_file: str, requested_mode: str | None = None):
    """Load node-loss counts from ``gd_loss_node_summary.tsv``.

    Returns
    -------
    tuple
        ``(counts, resolved_mode, used)`` where ``counts`` maps nodes to
        ``2_0/2_1/1_0/total`` counts.
    """
    requested_mode = normalize_visualizer_mode(requested_mode, default=None)
    with open(node_summary_file, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if not fieldnames:
        return {}, requested_mode or "parsimony", False

    explicit_modes = []
    for row in rows:
        mode_text = normalize_visualizer_mode(row.get("Mode"), default=None)
        if mode_text and mode_text not in explicit_modes:
            explicit_modes.append(mode_text)

    if explicit_modes:
        if requested_mode in explicit_modes:
            resolved_mode = requested_mode
        else:
            resolved_mode = explicit_modes[0]
            if requested_mode and requested_mode != resolved_mode:
                logger.info(
                    "Visualizer mode %s overridden by node summary mode %s from %s",
                    requested_mode,
                    resolved_mode,
                    node_summary_file,
                )
    elif {"C2_0", "C2_1", "C1_0"} <= set(fieldnames):
        resolved_mode = "accumulate"
    elif {"P2_0", "P2_1", "P1_0"} <= set(fieldnames):
        resolved_mode = "parsimony"
    else:
        resolved_mode = requested_mode or "parsimony"

    prefix = "C" if resolved_mode == "accumulate" else "P"
    c20_key = f"{prefix}2_0"
    c21_key = f"{prefix}2_1"
    c10_key = f"{prefix}1_0"

    counts = {}
    for row in rows:
        row_mode = normalize_visualizer_mode(row.get("Mode"), default=resolved_mode)
        if row_mode != resolved_mode:
            continue
        node_name = str(row.get("Node", "")).strip()
        if not node_name:
            continue
        c20 = int(row.get(c20_key, 0) or 0)
        c21 = int(row.get(c21_key, 0) or 0)
        c10 = int(row.get(c10_key, 0) or 0)
        total = int(row.get("Node loss count", c20 + c21 + c10) or 0)
        counts[node_name] = {
            "2_0": c20,
            "2_1": c21,
            "1_0": c10,
            "total": total,
        }
    return counts, resolved_mode, bool(counts)


def tree_has_numbered_internal_nodes(sptree) -> bool:
    """Return True when the species tree has S-numbered internal nodes."""
    return any(
        (not node.is_leaf()) and is_internal_node(str(node.name or ""))
        for node in sptree.traverse()
    )


def load_visualizer_species_tree(loss_result_file: str, species_tree_file: str):
    """Load the numbered species tree needed by GD_Loss_Visualizer.

    Tracker outputs use S-numbered internal nodes. If the user passes an
    unnumbered species tree, prefer the tracker-generated ``numed_sptree.nwk``
    from the same directory as ``gd_loss.csv`` so internal-node loss labels can
    be rendered.
    """
    from ete3 import Tree

    sptree = Tree(species_tree_file, format=1)
    if tree_has_numbered_internal_nodes(sptree):
        return sptree

    candidate = os.path.join(os.path.dirname(os.path.abspath(loss_result_file)), "numed_sptree.nwk")
    if os.path.exists(candidate):
        numbered_tree = Tree(candidate, format=1)
        if tree_has_numbered_internal_nodes(numbered_tree):
            logger.warning(
                "Input species tree has no S-numbered internal nodes; using tracker-generated numbered tree: %s",
                candidate,
            )
            return numbered_tree

    logger.warning(
        "Input species tree has no S-numbered internal nodes and no usable numed_sptree.nwk was found next to %s; "
        "internal-node GD-loss labels may be missing.",
        loss_result_file,
    )
    return sptree


def infer_visualizer_mode(filepath: str) -> str:
    """Infer how node-loss counts should be described in the figure.

    GD_Loss_Visualizer renders GD_Loss_Tracker output directly and does not
    ask the user to restate the counting mode. In practice, the tracker
    default is parsimony, so files without an explicit accumulate hint are
    interpreted as parsimony-style node counts.
    """
    node_summary_file = find_node_summary_file(filepath)
    if node_summary_file:
        _, resolved_mode, used = load_node_summary_counts(node_summary_file, requested_mode=None)
        if used:
            return resolved_mode
    if 'accumulate' in str(filepath).lower():
        return 'accumulate'
    return 'parsimony'


def get_event_stats(filepath: str):
    """Compute GD event counts and event-level 2_2/2_1/2_0 counts by GD node."""
    gd_birth_sets = defaultdict(set)
    best_event_type_by_level = defaultdict(dict)
    priority = {"2_2": 0, "2_1": 1, "2_0": 2}

    logger.info('Reading file: %s...', filepath)
    for cols, idx in _read_summary_rows(filepath):
        if len(cols) < 4:
            continue
        tree_i = _find_col(idx, "Tree ID", "tree_ID")
        gd_i = _find_col(idx, "GD ID", "gd_ID")
        level_i = _find_col(idx, "GD burst node", "GD Burst Node", "GD Event Node", "level")
        if tree_i is None or gd_i is None or level_i is None:
            continue
        tree_id = cols[tree_i]
        gd_id = cols[gd_i]
        level_node = cols[level_i]
        event_key = f"{tree_id}|{gd_id}"
        gd_birth_sets[level_node].add(event_key)
        loss_type_i = _find_col(idx, 'Species final loss type')
        if loss_type_i is None or loss_type_i >= len(cols):
            continue
        event_type = cols[loss_type_i].strip().replace("-", "_")
        if event_type not in priority:
            continue
        prev = best_event_type_by_level[level_node].get(event_key)
        if prev is None or priority[event_type] > priority[prev]:
            best_event_type_by_level[level_node][event_key] = event_type

    gd_births = {node: len(keys) for node, keys in gd_birth_sets.items()}
    event_loss_data = {}
    for node, event_map in best_event_type_by_level.items():
        counts = {'2_0': 0, '2_1': 0, '2_2': 0}
        for event_type in event_map.values():
            counts[event_type] += 1
        event_loss_data[node] = counts
    return gd_births, event_loss_data


def get_cumulative_counts_from_paths(filepath: str):
    """Aggregate cumulative path-level losses from ``loss_path``."""
    counts = defaultdict(lambda: {'2_0': 0, '2_1': 0, '1_0': 0, 'total': 0})
    for cols, idx in _read_summary_rows(filepath):
        path_i = _find_col(idx, 'GD loss path', 'GD Loss Path', 'Loss path', 'Loss Path', 'loss_path')
        if path_i is None or path_i >= len(cols):
            continue
        for node_name, loss_type in identify_loss_detail(cols[path_i].strip()):
            counts[node_name][loss_type] += 1
            counts[node_name]['total'] += 1
    return dict(counts)


def get_node_event_counts(filepath: str, mode: str = 'parsimony'):
    """Aggregate node-level loss events directly from tracker output.

    For parsimony-style summaries, tracker output may broadcast the same shared
    node event onto every descendant species row so the table stays readable.
    Visualization should collapse those back to one logical event per
    ``Tree ID`` + ``GD ID`` + ``node`` + ``type``. Accumulate summaries must not
    be collapsed because each row represents one counted species path.
    """
    mode = str(mode or 'parsimony').lower()
    counts = defaultdict(lambda: {'2_0': 0, '2_1': 0, '1_0': 0, 'total': 0})
    seen_event_keys = set()
    used = False
    for cols, idx in _read_summary_rows(filepath):
        tree_i = _find_col(idx, 'Tree ID', 'tree_ID')
        gd_i = _find_col(idx, 'GD ID', 'gd_ID')
        summary_i = _find_col(idx, 'Node loss events')
        if summary_i is None and mode == 'parsimony':
            summary_i = _find_col(idx, 'GD loss pattern')
        if summary_i is not None and summary_i < len(cols):
            parsed = parse_path_count_summary(cols[summary_i].strip())
            events = parse_node_loss_events(parsed["node_events"])
        else:
            continue
        if events:
            used = True
        tree_id = cols[tree_i] if tree_i is not None and tree_i < len(cols) else None
        gd_id = cols[gd_i] if gd_i is not None and gd_i < len(cols) else None
        for node_name, loss_type in events:
            if mode == 'parsimony' and tree_id is not None and gd_id is not None:
                event_key = (tree_id, gd_id, node_name, loss_type)
                if event_key in seen_event_keys:
                    continue
                seen_event_keys.add(event_key)
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
    """Fallback parsimony-like compression from species final loss assignments."""
    event_rows = {}
    for cols, idx in _read_summary_rows(filepath):
        loss_type_i = _find_col(idx, 'Species final loss type')
        species_i = _find_col(idx, 'Species', 'species')
        level_i = _find_col(idx, 'GD burst node', 'GD Burst Node', 'GD Event Node', 'level')
        if (
            loss_type_i is None
            or species_i is None
            or level_i is None
            or max(loss_type_i, species_i, level_i) >= len(cols)
        ):
            continue
        tree_i = _find_col(idx, 'Tree ID', 'tree_ID')
        gd_i = _find_col(idx, 'GD ID', 'gd_ID')
        if tree_i is None or gd_i is None:
            continue
        event_key = (cols[tree_i], cols[gd_i])
        slot = event_rows.setdefault(event_key, {
            'level': cols[level_i],
            'loss_one_species': set(),
            'loss_two_species': set(),
        })
        loss_type = cols[loss_type_i].strip().replace("-", "_")
        species = cols[species_i]
        if loss_type == '2_0':
            slot['loss_two_species'].add(species)
        elif loss_type == '2_1':
            slot['loss_one_species'].add(species)

    parsimony_counts = defaultdict(lambda: {'2_0': 0, '2_1': 0, '1_0': 0, 'total': 0})
    for event in event_rows.values():
        try:
            gd_node = sptree & event['level']
        except Exception:
            continue
        loss_two_cover = _get_maximal_cover_nodes(gd_node, event['loss_two_species'])
        covered_by_two = set()
        for node_name in loss_two_cover:
            parsimony_counts[node_name]['1_0'] += 1
            parsimony_counts[node_name]['total'] += 1
            try:
                covered_by_two.update((sptree & node_name).get_leaf_names())
            except Exception:
                pass
        residual_loss_one = set(event['loss_one_species']) - covered_by_two
        for node_name in _get_maximal_cover_nodes(gd_node, residual_loss_one):
            parsimony_counts[node_name]['2_1'] += 1
            parsimony_counts[node_name]['total'] += 1
    return dict(parsimony_counts)


def format_event_birth_summary(gd_birth_count: int) -> str:
    return f"B={gd_birth_count}"


def format_node_loss_breakdown(node_stats: dict, loss_label: str) -> str:
    """Format per-node loss-type breakdown for display on the tree."""
    return (
        f"{loss_label}20/21/10="
        f"{node_stats.get('2_0', 0)}/"
        f"{node_stats.get('2_1', 0)}/"
        f"{node_stats.get('1_0', 0)}"
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
        from ete3 import TreeStyle, NodeStyle, TextFace
    except ImportError:
        logger.error('ete3 is required for visualization. Install with: pip install ete3')
        return

    visualizer_mode = infer_visualizer_mode(filepath)
    gd_births, event_loss_data = get_event_stats(filepath)

    node_loss_source = filepath
    node_loss_counts = {}
    node_summary_file = find_node_summary_file(filepath)
    if node_summary_file:
        node_loss_counts, visualizer_mode, used_node_summary = load_node_summary_counts(
            node_summary_file,
            requested_mode=visualizer_mode,
        )
        if used_node_summary:
            node_loss_source = node_summary_file
    else:
        used_node_summary = False

    if not used_node_summary:
        node_loss_counts, has_node_events = get_node_event_counts(filepath, visualizer_mode)
        if has_node_events:
            node_loss_source = filepath
    else:
        has_node_events = True

    if not node_loss_counts:
        if visualizer_mode == 'accumulate':
            node_loss_counts = get_cumulative_counts_from_paths(filepath)
            node_loss_source = f"{filepath} [GD loss path fallback]"
        else:
            node_loss_counts = get_fallback_parsimony_counts(filepath, sptree)
            node_loss_source = f"{filepath} [species-loss fallback]"

    if visualizer_mode == 'accumulate':
        loss_label = 'C'
        loss_desc = 'accumulate node-loss count from GD_Loss_Tracker'
        face_color = 'darkorange'
    else:
        loss_label = 'P'
        loss_desc = 'parsimony node-loss count from GD_Loss_Tracker (shared MRCA losses counted once)'
        face_color = 'darkgreen'
    logger.info('Node-loss source: %s', node_loss_source)
    logger.info('Visualizer mode: %s', visualizer_mode)

    sptree.ladderize()
    sptree.sort_descendants('support')

    if visualizer_mode == 'accumulate':
        logger.info(
            '\n%-15s | %-8s | %-10s | %-8s | %-8s | %-8s',
            'Node', 'GDcount', 'NodeLoss', 'C2_0', 'C2_1', 'C1_0'
        )
        logger.info('-' * 76)
    else:
        logger.info(
            '\n%-15s | %-8s | %-8s | %-8s | %-10s | %-10s',
            'Node', 'GDcount', 'E2_0', 'E2_1', 'NodeLoss', 'Label'
        )
        logger.info('-' * 72)

    for node in sptree.traverse():
        nstyle = NodeStyle()
        nstyle['size'] = 0
        nstyle['hz_line_width'] = 1
        nstyle['vt_line_width'] = 1
        nstyle['hz_line_type'] = 0
        nstyle['vt_line_type'] = 0
        nstyle['hz_line_color'] = 'black'
        nstyle['vt_line_color'] = 'black'
        node.set_style(nstyle)
        node_name = node.name
        is_numbered_internal = (not node.is_leaf()) and is_internal_node(str(node_name or ""))

        if node_name in gd_births or node_name in event_loss_data or node_name in node_loss_counts:
            if visualizer_mode == 'accumulate' and not is_numbered_internal:
                continue
            event_stats = event_loss_data.get(node_name, {'2_0': 0, '2_1': 0, '2_2': 0})
            node_stats = node_loss_counts.get(node_name, {'2_0': 0, '2_1': 0, '1_0': 0, 'total': 0})
            gd_birth_count = gd_births.get(node_name, 0)
            node_total = node_stats.get('total', 0)
            if is_numbered_internal and gd_birth_count > 0:
                node.add_face(
                    TextFace(f"{node_name} {format_event_birth_summary(gd_birth_count)}", fsize=6, fgcolor='blue'),
                    column=0,
                    position='branch-top',
                )

            if node_total > 0 and not node.is_leaf():
                node.add_face(TextFace(f"{loss_label}={node_total}", fsize=6, fgcolor=face_color), column=0, position='branch-bottom')
                node.add_face(
                    TextFace(format_node_loss_breakdown(node_stats, loss_label), fsize=6, fgcolor=face_color),
                    column=0,
                    position='branch-bottom',
                )

            if visualizer_mode == 'accumulate':
                if gd_birth_count > 0 or node_total > 0:
                    logger.info(
                        '%-15s | %-8d | %-10d | %-8d | %-8d | %-8d',
                        node_name,
                        gd_birth_count,
                        node_total,
                        node_stats.get('2_0', 0),
                        node_stats.get('2_1', 0),
                        node_stats.get('1_0', 0),
                    )
            elif node_total > 0 or gd_birth_count > 0:
                logger.info(
                    '%-15s | %-8d | %-8d | %-8d | %-10d | %-10s',
                    node_name,
                    gd_birth_count,
                    event_stats.get('2_0', 0),
                    event_stats.get('2_1', 0),
                    node_total,
                    loss_label,
                )

    ts = TreeStyle()
    ts.scale = 40
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.show_branch_length = False
    ts.show_border = False
    ts.extra_branch_line_type = 0
    ts.extra_branch_line_color = 'black'
    ts.complete_branch_lines_when_necessary = True
    ts.draw_guiding_lines = True
    ts.guiding_lines_type = 0
    ts.guiding_lines_color = 'black'
    ts.branch_vertical_margin = 3
    ts.margin_right = 20
    ts.legend_position = 1
    ts.title.add_face(TextFace('Legend:', bold=True, fsize=10), column=0)
    ts.title.add_face(TextFace('  Sxx=node identifier; B=GD burst count', fsize=6, fgcolor='blue'), column=0)
    ts.title.add_face(TextFace(f'  {loss_label}: {loss_desc}', fsize=6, fgcolor=face_color), column=0)
    ts.title.add_face(TextFace(f'  {loss_label}20/21/10: node-loss breakdown by 2_0, 2_1, and 1_0 counts', fsize=6, fgcolor=face_color), column=0)

    try:
        sptree.convert_to_ultrametric()
    except Exception:
        pass

    for leaf in sptree.iter_leaves():
        leaf.img_style['hz_line_width'] = 1
        leaf.img_style['vt_line_width'] = 1
        leaf.img_style['hz_line_type'] = 0
        leaf.img_style['vt_line_type'] = 0
        leaf.img_style['hz_line_color'] = 'black'
        leaf.img_style['vt_line_color'] = 'black'
        leaf_stats = node_loss_counts.get(leaf.name, {'2_0': 0, '2_1': 0, '1_0': 0, 'total': 0})
        leaf_total = leaf_stats.get('total', 0)
        if leaf_total > 0:
            total_face = TextFace(f"{loss_label}={leaf_total}", fsize=6, fgcolor=face_color)
            leaf.add_face(total_face, column=0, position='branch-bottom')

            breakdown_face = TextFace(
                format_node_loss_breakdown(leaf_stats, loss_label),
                fsize=6,
                fgcolor=face_color,
            )
            leaf.add_face(breakdown_face, column=0, position='branch-bottom')

        name_face = TextFace(leaf.name, fsize=7, fgcolor='black')
        name_face.margin_left = 3
        leaf.add_face(name_face, column=2, position='aligned')

    sptree.render(output_file, w=260, units='mm', tree_style=ts)
    logger.info('Visualization saved to %s', output_file)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Visualize gene duplication losses on a species tree.')
    parser.add_argument('loss_summary_file', help='Path to one GD loss summary file.')
    parser.add_argument('species_tree_file', help='Path to the species tree file (Newick format).')
    parser.add_argument('-o', '--output', default='gd_loss_pie_visualizer.PDF', help='Output PDF file path (default: gd_loss_pie_visualizer.PDF).')
    args = parser.parse_args()

    sptree = load_visualizer_species_tree(args.loss_summary_file, args.species_tree_file)
    visualizer_sptree(
        args.loss_summary_file,
        sptree,
        output_file=args.output,
    )
