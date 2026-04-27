"""
Gene duplication loss tracking for the PhyloTracer pipeline.

This module evaluates loss events following duplication nodes, summarizes
loss paths, and generates tabular reports for downstream analyses.
"""

import logging
import re
from typing import Optional

import os
from collections import defaultdict

logger = logging.getLogger(__name__)

import pandas as pd
from ete3 import PhyloTree

from phylotracer import (
    rename_input_tre,
    get_species_set,
    num_sptree,
    read_and_return_dict,
    find_dup_node,
    num_tre_node,
    annotate_gene_tree,
    gene_id_transfer,
    map_species_set_to_node,
    parse_loss_path,
    is_internal_node,
    parse_node_number,
)

# ======================================================
# Constants
# ======================================================

LOSS_TYPE_2_2 = "2_2"  # Both copies retained
LOSS_TYPE_2_1 = "2_1"  # One copy lost
LOSS_TYPE_2_0 = "2_0"  # Absent from both GD child subclades
LOSS_TYPE_1_0 = "1_0"  # Single copy lost
LOSS_TYPE_NA = "NA"

SORT_SENTINEL = 10**9  # Sentinel for non-S-node sorting

# ======================================================
# Helper Functions
# ======================================================


def node_sort_key(name: str) -> int:
    """Return sort key for internal node names (S0, S1, etc)."""
    num = parse_node_number(name)
    return num if num is not None else SORT_SENTINEL


def classify_copy_transition(prev_copy: int, curr_copy: int) -> str:
    """Classify copy number transition into loss type."""
    if curr_copy >= prev_copy:
        return LOSS_TYPE_NA
    if prev_copy == 2 and curr_copy == 0:
        return LOSS_TYPE_2_0
    elif prev_copy == 2 and curr_copy == 1:
        return LOSS_TYPE_2_1
    elif prev_copy == 1 and curr_copy == 0:
        return LOSS_TYPE_1_0
    else:
        return LOSS_TYPE_2_0 if (prev_copy - curr_copy) >= 2 else LOSS_TYPE_2_1


# ======================================================
# Section 1: Duplication Loss Validation
# ======================================================


def legacy_is_valid_duplication_loss(species_voucher, dup_node, renamed_sptree):
    """
    Deprecated legacy validator for duplication loss.

    Parameters
    ----------
    species_voucher : object
        Voucher identifier for the target species.
    dup_node : object
        Duplication node in the gene tree.
    renamed_sptree : object
        Species tree with renamed labels for mapping.

    Returns
    -------
    bool
        Legacy boolean flag; retained only for backward compatibility.

    Assumptions
    -----------
    Deprecated: this function uses a strict AND criterion (left and right
    both present before calling loss), which is inconsistent with event-level
    2_2/2_1/2_0 classification and should not be used in new code.
    """
    gd_species_set = get_species_set(dup_node)

    if len(gd_species_set) == 1:
        mapped_node = renamed_sptree & list(gd_species_set)[0]
    else:
        mapped_node = renamed_sptree.get_common_ancestor(gd_species_set)

    if mapped_node is None:
        return False

    if species_voucher not in set(mapped_node.get_leaf_names()):
        return False

    children = dup_node.get_children()
    if len(children) != 2:
        return False

    left, right = children

    left_species = {leaf.split("_")[0] for leaf in left.get_leaf_names()}
    right_species = {leaf.split("_")[0] for leaf in right.get_leaf_names()}

    should_have_two = (species_voucher in left_species) and (species_voucher in right_species)

    if not should_have_two:
        return False

    observed_copy_num = 0
    for leaf in dup_node.get_leaf_names():
        if leaf.split("_")[0] == species_voucher:
            observed_copy_num += 1

    return observed_copy_num < 2


def classify_species_copy_state(
    species_voucher: str,
    dup_node: object,
    mapped_node: object,
    family_species_set: Optional[set] = None,
    include_unobserved_species: bool = False,
) -> tuple:
    """Classify post-GD copy state for one species at one duplication node.

    Parameters
    ----------
    species_voucher : str
        Species voucher to classify.
    dup_node : object
        Duplication node with two child lineages.
    mapped_node : object
        Species-tree node mapped from the duplication event.
    family_species_set : set, optional
        Species vouchers observed in the current gene family.

    Returns
    -------
    tuple
        ``(loss_type, left_has, right_has, confidence)`` where ``loss_type`` is
        one of ``2_2``, ``2_1``, ``2_0``, or ``NA``.
    """
    children = dup_node.get_children()
    if mapped_node is None or len(children) != 2:
        return LOSS_TYPE_NA, False, False, "invalid_node"

    if species_voucher not in set(mapped_node.get_leaf_names()):
        return LOSS_TYPE_NA, False, False, "out_of_mapped_node"

    left_child, right_child = children
    left_species = get_species_set(left_child)
    right_species = get_species_set(right_child)

    left_has = species_voucher in left_species
    right_has = species_voucher in right_species

    if family_species_set is not None and species_voucher not in family_species_set:
        if not include_unobserved_species:
            return LOSS_TYPE_NA, left_has, right_has, "unobserved_in_family"
        confidence = "unobserved_in_family_included"
    else:
        confidence = "observable"

    if left_has and right_has:
        loss_type = LOSS_TYPE_2_2
    elif left_has or right_has:
        loss_type = LOSS_TYPE_2_1
    else:
        loss_type = LOSS_TYPE_2_0

    return loss_type, left_has, right_has, confidence


def summarize_small_loss_types_from_path(path_str: str) -> tuple:
    """Extract small loss-type transitions from a path string.

    Parameters
    ----------
    path_str : str
        Path string formatted as ``node(count)->node(count)``.

    Returns
    -------
    tuple
        ``(path_count_types, c20, c21, c10, transition_keys)`` where counts are path-count
        transitions (not biological copy-number states).
    """
    parsed = parse_loss_path(path_str)
    if not parsed or len(parsed) < 2:
        return LOSS_TYPE_NA, 0, 0, 0, []

    c20 = c21 = c10 = 0
    ordered_types = []
    transition_keys = []
    for i in range(1, len(parsed)):
        prev_copy = parsed[i - 1][1]
        curr_node, curr_copy = parsed[i]

        loss_type = classify_copy_transition(prev_copy, curr_copy)
        if loss_type == LOSS_TYPE_NA:
            continue

        ordered_types.append(loss_type)
        transition_keys.append((curr_node, loss_type))

        if loss_type == LOSS_TYPE_2_0:
            c20 += 1
        elif loss_type == LOSS_TYPE_2_1:
            c21 += 1
        elif loss_type == LOSS_TYPE_1_0:
            c10 += 1

    return (",".join(ordered_types) if ordered_types else LOSS_TYPE_NA), c20, c21, c10, transition_keys


def build_transition_summary(transition_keys):
    """Convert transition-key list into exported summary fields."""
    if not transition_keys:
        return LOSS_TYPE_NA, 0, 0, 0, LOSS_TYPE_NA

    ordered_types = [trans_type for _, trans_type in transition_keys]
    c20 = ordered_types.count(LOSS_TYPE_2_0)
    c21 = ordered_types.count(LOSS_TYPE_2_1)
    c10 = ordered_types.count(LOSS_TYPE_1_0)
    node_events = ";".join(f"{node_name}|{trans_type}" for node_name, trans_type in transition_keys)
    return ",".join(ordered_types), c20, c21, c10, node_events


def build_path_count_summary(node_events, types_str, c20, c21, c10):
    """Backward-compatible helper retained for older callers."""
    if not types_str or types_str == LOSS_TYPE_NA:
        return LOSS_TYPE_NA
    return str(types_str).replace(",", ";")


def build_path_count_summary_by_mode(mode, transition_keys):
    """Build one concise node-event summary from mode-specific transition keys.

    ``transition_keys`` has already been prepared according to the selected mode:
    accumulate keeps every path transition, while parsimony keeps collapsed shared
    transitions. Export the same node:type representation for both modes so users
    can verify node counts from one column.
    """
    if not transition_keys:
        return LOSS_TYPE_NA
    return ";".join(f"{node_name}:{loss_type}" for node_name, loss_type in transition_keys)


def get_path_node_order(path_str: str) -> dict:
    """Return node order on one loss path for stable event display."""
    parsed = parse_loss_path(path_str)
    return {node_name: idx for idx, (node_name, _) in enumerate(parsed)}


def infer_loss_node_from_path(path_str: str) -> str:
    """Return the first node where copy number decreases on a loss path."""
    parsed = parse_loss_path(path_str)
    if len(parsed) < 2:
        return LOSS_TYPE_NA

    prev = parsed[0][1]
    for node_name, curr in parsed[1:]:
        if curr < prev:
            return node_name
        prev = curr
    return LOSS_TYPE_NA


def identify_loss_detail(path_str: str):
    """Parse one loss path and return node-level copy-loss transitions.

    Returns tuples of ``(node_name, loss_type)`` where ``loss_type`` is one of
    ``2_0``, ``2_1`` or ``1_0``. The node is the node reached after the copy
    number decreases.
    """
    parsed_steps = parse_loss_path(path_str)
    if not parsed_steps:
        return []

    losses = []
    for idx in range(1, len(parsed_steps)):
        _, prev_copy = parsed_steps[idx - 1]
        curr_node, curr_copy = parsed_steps[idx]

        loss_type = classify_copy_transition(prev_copy, curr_copy)
        if loss_type != LOSS_TYPE_NA:
            losses.append((curr_node, loss_type))
    return losses


def infer_parsimony_transition_keys(
    rows,
    mapped_clade,
    voucher2taxa_dic,
    min_support_ratio: float = 0.5,
    min_support_species: int = 2,
):
    """Infer minimal ancestor-level loss events for one GD event.

    The rule is clade-cover based: if all descendant species of a species-tree node
    share a loss transition, record the loss once at that ancestor rather than once
    for every terminal species.
    """
    transition_types = (LOSS_TYPE_2_0, LOSS_TYPE_2_1, LOSS_TYPE_1_0)
    assigned = {row["species_voucher"]: [] for row in rows}
    def cover_nodes(node, affected_species, out_nodes):
        leaves = set(node.get_leaf_names())
        covered = leaves & affected_species
        if not covered:
            return
        if leaves <= affected_species:
            out_nodes.append(node)
            return
        children = node.get_children()
        if not children:
            out_nodes.append(node)
            return
        for child in children:
            cover_nodes(child, affected_species, out_nodes)

    row_by_species = {row["species_voucher"]: row for row in rows}
    for trans_type in transition_types:
        affected_species = {
            row["species_voucher"]
            for row in rows
            if any(t == trans_type for _, t in row["raw_transition_keys"])
        }
        if not affected_species:
            continue
        selected_nodes = []
        cover_nodes(mapped_clade, affected_species, selected_nodes)
        for node in selected_nodes:
            desc = set(node.get_leaf_names()) & affected_species
            if not desc:
                continue
            support_species = len(desc)
            node_leaves = set(node.get_leaf_names())
            local_species_total = len(node_leaves)
            support_ratio = (support_species / local_species_total) if local_species_total else 0.0
            if not node.is_leaf():
                if support_species < min_support_species:
                    continue
                if support_ratio < min_support_ratio:
                    continue
            node_label = voucher2taxa_dic.get(node.name, node.name)
            for species_voucher in desc:
                assigned[species_voucher].append((node_label, trans_type))

    for row in rows:
        node_order = row["path_node_order"]
        row["effective_transition_keys"] = sorted(
            assigned.get(row["species_voucher"], []),
            key=lambda item: (node_order.get(item[0], SORT_SENTINEL), item[0], item[1]),
        )

    return rows


# ======================================================
# Section 2: Mapping and Path Utilities
# ======================================================


def get_maptree_internal_node_name_set(node, sptree):
    sps = get_species_set(node)

    node1 = node.up
    map_nodename1 = get_maptree_name(sptree, get_species_set(node1))
    map_node1 = sptree & map_nodename1

    names = []
    for i in sps:
        clade = sptree & i
        path = get_two_nodes_path_str(clade, map_node1)
        names += path

    # Node copy-state on path annotations must stay in biological range {0,1,2},
    # so each child-side clade contributes at most one vote per internal node.
    return set(names)


def get_two_subclade_maptree_node_name_lst(max_clade, sptree):
    clade_up = max_clade.get_children()[0]
    clade_up_set = get_maptree_internal_node_name_set(clade_up, sptree)
    clade_down = max_clade.get_children()[1]
    clade_down_set = get_maptree_internal_node_name_set(clade_down, sptree)
    up_down_lst = list(clade_up_set) + list(clade_down_set)
    return up_down_lst


def get_maptree_internal_node_name_count_dic(
    max_clade,
    max_clade2sp,
    sptree,
):
    up_down_lst = get_two_subclade_maptree_node_name_lst(max_clade, sptree)
    dic = {i.name: 0 for i in max_clade2sp.traverse()}
    for i in up_down_lst:
        if i in dic:
            dic[i] += 1
    keys_with_zero_value = [
        key
        for key, value in dic.items()
        if value == 0 and key not in sptree.get_leaf_names()
    ]
    return dic, keys_with_zero_value


def get_maptree_name(sptree, sps_set):
    if len(sps_set) != 1:
        com = sptree.get_common_ancestor(sps_set)
        return com.name
    else:
        com = sptree & list(sps_set)[0]
        return com.name


def get_two_nodes_path_str(start_node, end_node):
    nodes = []
    current_node = start_node
    while current_node != end_node:
        nodes.append(current_node.name)
        current_node = current_node.up
        if current_node is None:
            break
    nodes.append(end_node.name)
    re_nodes = list(reversed(nodes))

    return re_nodes


def get_tips_to_clade_path_lst(taget_node: object, dic) -> list:
    path_str_lst = []
    for i in taget_node:
        path_str = get_two_nodes_path_str(i, taget_node)
        numlist = [dic[k] for k in path_str]
        result = "->".join([f"{name}({count})" for name, count in zip(path_str, numlist)])

        path_str_lst.append(result)
    return path_str_lst


def get_maptree_node_count_dic(sp_list, map_clade):
    sp_dic = {i.name: 0 for i in map_clade.traverse()}
    for i in sp_list:
        if i in sp_dic:
            sp_dic[i] += 1

    tips = set([j for j in sp_dic.keys() if not (j.startswith("S") and j[1:].isdigit())])  # only match S+digits as internal nodes, not species like Solanum
    for k, v in sp_dic.items():
        if k.startswith("S") and k[1:].isdigit():  # match S+digits only
            if v == 0:
                t = map_clade & k
                sp = get_species_set(t)
                if len(sp & tips) != 0:
                    sp_dic[k] = 1
            elif v == 1:
                t = map_clade & k
                sp = get_species_set(t)
                if len(sp & tips) != 0:
                    jiao = list(sp & tips)[0]
                    if sp_dic[jiao] == 2:
                        sp_dic[k] = 2
    return sp_dic


# ======================================================
# Section 3: Gene Pair and Path Extraction
# ======================================================


def get_path_str_with_count_num_lst(
    tre_id,
    gd_id_start,
    genetree,
    renamed_sptree,
    new_named_gene2gene_dic,
    voucher2taxa_dic,
    include_unobserved_species=False,
    node_count_mode: str = "parsimony",
    parsimony_min_support_ratio: float = 0.5,
    parsimony_min_support_species: int = 2,
):
    """
    Process duplication nodes and record gene pairs and loss paths.

    Parameters
    ----------
    tre_id : object
        Tree identifier.
    gd_id_start : int
        Starting GD identifier for numbering.
    genetree : object
        Gene tree containing duplication nodes.
    renamed_sptree : object
        Renamed species tree used for mapping.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene identifiers.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.

    Returns
    -------
    tuple
        (records, path_info_list, next_gd_id).

    Assumptions
    -----------
    Duplication nodes can be detected by ``find_dup_node``.
    """
    records = []
    path_info_list = []
    dup_nodes = find_dup_node(genetree, renamed_sptree)

    gd_id = gd_id_start

    for dup_node in dup_nodes:
        sp = get_species_set(dup_node)

        max_clade2sp = map_species_set_to_node(renamed_sptree,sp)
        
        gd_node_name = max_clade2sp.name

        voucher_to_pretty_path = {}
        if len(sp) > 1:
            count_dic, _ = get_maptree_internal_node_name_count_dic(
                dup_node,
                max_clade2sp,
                renamed_sptree,
            )
            path_str_lst = get_tips_to_clade_path_lst(max_clade2sp, count_dic)

            for path_str in path_str_lst:
                path_info_list.append((path_str, gd_id, gd_node_name))

            voucher_to_raw_path = {}
            for path_str in path_str_lst:
                last_voucher = path_str.split("->")[-1].split("(")[0].strip()
                voucher_to_raw_path[last_voucher] = path_str

            for voucher, raw_path in voucher_to_raw_path.items():
                parts = raw_path.split("->")
                last_part = parts[-1]
                count_part = last_part.partition("(")[2]
                taxa_name = voucher2taxa_dic.get(voucher, voucher)
                parts[-1] = (taxa_name + "(" + count_part) if count_part else taxa_name
                voucher_to_pretty_path[voucher] = "->".join(parts)

        gd_level_name = voucher2taxa_dic.get(max_clade2sp.name, max_clade2sp.name)
        family_species_set = get_species_set(genetree)
        children = dup_node.get_children()
        if len(children) != 2:
            gd_id += 1
            continue

        left_child, right_child = children
        left_map = {}
        right_map = {}

        for gene_name in left_child.get_leaf_names():
            species = gene_name.split("_", 1)[0]
            left_map.setdefault(species, []).append(gene_name)
        for gene_name in right_child.get_leaf_names():
            species = gene_name.split("_", 1)[0]
            right_map.setdefault(species, []).append(gene_name)

        clade_species = set(max_clade2sp.get_leaf_names())
        seen_species = set()
        seen_keys = set()
        gd_counter = {LOSS_TYPE_2_2: 0, LOSS_TYPE_2_1: 0, LOSS_TYPE_2_0: 0, "missing_data": 0}
        mode_normalized = node_count_mode
        pending_rows = []
        for species_voucher in sorted(clade_species):
            if species_voucher in seen_species:
                continue
            seen_species.add(species_voucher)
            rec_key = (tre_id, gd_id, species_voucher)
            if rec_key in seen_keys:
                continue
            seen_keys.add(rec_key)

            loss_type, left_has, right_has, confidence = classify_species_copy_state(
                species_voucher,
                dup_node,
                max_clade2sp,
                family_species_set,
                include_unobserved_species=include_unobserved_species,
            )
            left_genes = sorted(left_map.get(species_voucher, []))
            right_genes = sorted(right_map.get(species_voucher, []))
            orig_left_genes = [new_named_gene2gene_dic.get(g, g) for g in left_genes]
            orig_right_genes = [new_named_gene2gene_dic.get(g, g) for g in right_genes]

            orig_a = ",".join(orig_left_genes) if orig_left_genes else LOSS_TYPE_NA
            orig_b = ",".join(orig_right_genes) if orig_right_genes else LOSS_TYPE_NA
            pretty_path = voucher_to_pretty_path.get(species_voucher, LOSS_TYPE_NA)
            _, _, _, _, transition_keys = summarize_small_loss_types_from_path(pretty_path)

            if confidence == "unobserved_in_family" and not include_unobserved_species:
                major_loss_class = "missing_data"
            elif loss_type == LOSS_TYPE_2_0:
                major_loss_class = "candidate_absent_both_subclades"
            elif loss_type == LOSS_TYPE_2_1:
                major_loss_class = "loss_one_copy"
            elif loss_type == LOSS_TYPE_2_2:
                major_loss_class = "no_loss"
            else:
                major_loss_class = LOSS_TYPE_NA

            if loss_type in gd_counter:
                gd_counter[loss_type] += 1
            elif major_loss_class == "missing_data":
                gd_counter["missing_data"] += 1

            pending_rows.append(
                {
                    "tre_id": tre_id,
                    "gd_id": gd_id,
                    "support": dup_node.support,
                    "gd_level_name": gd_level_name,
                    "species_display": voucher2taxa_dic.get(species_voucher, species_voucher),
                    "species_voucher": species_voucher,
                    "gene1": orig_a,
                    "gene2": orig_b,
                    "left_gene_n": len(left_genes),
                    "right_gene_n": len(right_genes),
                    "loss_path": pretty_path,
                    "gd_node_name": gd_node_name,
                    "loss_type": loss_type,
                    "left_has": int(left_has),
                    "right_has": int(right_has),
                    "loss_confidence": confidence,
                    "major_loss_class": major_loss_class,
                    "raw_transition_keys": list(transition_keys),
                    "path_node_order": get_path_node_order(pretty_path),
                }
            )

        if mode_normalized == "parsimony":
            pending_rows = infer_parsimony_transition_keys(
                pending_rows,
                max_clade2sp,
                voucher2taxa_dic,
                min_support_ratio=parsimony_min_support_ratio,
                min_support_species=parsimony_min_support_species,
            )
        else:
            for row in pending_rows:
                row["effective_transition_keys"] = row["raw_transition_keys"]

        for row in pending_rows:
            effective_keys = row.get("effective_transition_keys", [])
            row["effective_transition_keys_for_export"] = effective_keys
            path_count_types, path_count_2_0, path_count_2_1, path_count_1_0, path_count_node_events = build_transition_summary(
                effective_keys
            )
            row["path_count_node_events"] = path_count_node_events
            row["path_count_types"] = path_count_types
            row["path_count_2_0"] = path_count_2_0
            row["path_count_2_1"] = path_count_2_1
            row["path_count_1_0"] = path_count_1_0
            row.pop("raw_transition_keys", None)
            row.pop("path_node_order", None)
            row.pop("effective_transition_keys", None)
            records.append(row)

        clade_species_count = len(clade_species)
        observed_species_in_clade = len(clade_species & family_species_set)
        n22 = gd_counter[LOSS_TYPE_2_2]
        n21 = gd_counter[LOSS_TYPE_2_1]
        n20 = gd_counter[LOSS_TYPE_2_0]
        n_missing = gd_counter["missing_data"]
        if (n22 + n21 + n20 + n_missing) != clade_species_count:
            logger.warning(
                "Species accounting mismatch: tre=%s, gd=%s, "
                "n22=%d, n21=%d, n20=%d, n_missing=%d, "
                "clade_species_count=%d",
                tre_id, gd_id, n22, n21, n20, n_missing, clade_species_count
            )
        if (
            not include_unobserved_species
            and (n22 + n21 + n20) != observed_species_in_clade
        ):
            logger.warning(
                "Observed-species accounting mismatch: tre=%s, gd=%s, "
                "n22=%d, n21=%d, n20=%d, "
                "observed_species_in_clade=%d",
                tre_id, gd_id, n22, n21, n20, observed_species_in_clade
            )

        gd_id += 1

    return records, path_info_list, gd_id


# ======================================================
# Section 4: Species Tree Numbering and Path Grouping
# ======================================================


def split_dict_by_first_last_char(original_dict):
    split_dicts = {}

    for key, value in original_dict.items():
        first_char = key.split("->")[0].split("(")[0]
        last_char = key.split("->")[-1].split("(")[0]

        new_key = f"{first_char}_{last_char}"

        if new_key not in split_dicts:
            split_dicts[new_key] = {}

        split_dicts[new_key][key] = value

    return split_dicts


# ======================================================
# Section 5: Loss Path Summaries
# ======================================================


def get_path_str_num_dic(
    tre_dic,
    sptree,
    gene2new_named_gene_dic,
    new_named_gene2gene_dic,
    voucher2taxa_dic,
    taxa2voucher_dic,
    target_species_list=None,
    allowed_gd_species_sets=None,
    include_unobserved_species=False,
    node_count_mode: str = "parsimony",
    parsimony_min_support_ratio: float = 0.5,
    parsimony_min_support_species: int = 2,
):
    """
    Aggregate loss-path statistics with optional filtering constraints.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to file paths.
    sptree : object
        Species tree used for mapping.
    gene2new_named_gene_dic : dict
        Mapping from gene identifiers to renamed identifiers.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene identifiers.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.
    taxa2voucher_dic : dict
        Mapping from taxa labels to voucher identifiers.
    target_species_list : list, optional
        List of taxa to restrict loss endpoints.
    allowed_gd_species_sets : set, optional
        Set of frozensets representing allowed MRCA species sets.

    Returns
    -------
    tuple
        (path_str_num_dic, path2_treeid_dic).

    Assumptions
    -----------
    Species tree and gene tree identifiers are mutually compatible.
    """
    renamed_sptree = rename_input_tre(sptree, taxa2voucher_dic)
    path2_treeid_dic = {}
    path_str_num_dic = {}

    all_records = []
    all_path_info_list = []
    global_gd_id = 1

    for tre_id, tre_path in tre_dic.items():
        t = PhyloTree(tre_path)
        
        if len(t.children) != 2:
            logger.warning("%s is not a binary tree, skipping.", tre_id)
            continue
        t1 = rename_input_tre(t, gene2new_named_gene_dic)
        annotate_gene_tree(t1,renamed_sptree)
        records, path_info_list, new_gd_id = get_path_str_with_count_num_lst(
            tre_id,
            global_gd_id,
            t1,
            renamed_sptree,
            new_named_gene2gene_dic,
            voucher2taxa_dic,
            include_unobserved_species=include_unobserved_species,
            node_count_mode=node_count_mode,
            parsimony_min_support_ratio=parsimony_min_support_ratio,
            parsimony_min_support_species=parsimony_min_support_species,
        )
        global_gd_id = new_gd_id

        all_records.extend(records)
        for p_str, g_id, g_node in path_info_list:
            all_path_info_list.append((p_str, g_id, g_node, tre_id))

    target_voucher_set = set()
    if target_species_list:
        for taxa in target_species_list:
            vouchers = [k for k, v in voucher2taxa_dic.items() if v == taxa]
            target_voucher_set.update(vouchers)
        logger.info(
            "[Filter] Loss endpoint restricted to: %s (vouchers: %s)",
            target_species_list, sorted(target_voucher_set)
        )

    if allowed_gd_species_sets:
        logger.info(
            "[Filter] GD event must occur at EXACT MRCA node(s), covering "
            "%d specified clade(s).", len(allowed_gd_species_sets)
        )
    else:
        logger.info("[Filter] No restriction on GD node location.")

    filtered_records = []
    for rec in all_records:
        if target_species_list and rec["species_voucher"] not in target_voucher_set:
            continue

        if allowed_gd_species_sets:
            gd_node = renamed_sptree & rec["gd_node_name"]
            gd_leaves_voucher = gd_node.get_leaf_names()
            gd_leaves_taxa = frozenset(
                voucher2taxa_dic.get(leaf, leaf) for leaf in gd_leaves_voucher
            )
            if gd_leaves_taxa not in allowed_gd_species_sets:
                continue

        filtered_records.append(rec)

    for path_str, gd_id, gd_node_name, tre_id in all_path_info_list:
        last_node = path_str.split("->")[-1].split("(")[0].strip()

        if target_species_list:
            if last_node not in target_voucher_set:
                continue

        if allowed_gd_species_sets:
            gd_node = renamed_sptree & gd_node_name
            gd_leaves_voucher = gd_node.get_leaf_names()
            gd_leaves_taxa = frozenset(
                voucher2taxa_dic.get(leaf, leaf) for leaf in gd_leaves_voucher
            )
            if gd_leaves_taxa not in allowed_gd_species_sets:
                continue

        path_str_num_dic[path_str] = path_str_num_dic.get(path_str, 0) + 1
        path2_treeid_dic.setdefault(path_str, []).append(f"{tre_id}-{gd_id}")

    if not path_str_num_dic:
        logger.warning("No paths matched the filters. Result is empty.")

    return path_str_num_dic, path2_treeid_dic, filtered_records


# ======================================================
# Section 6: Sorting and Excel Reporting
# ======================================================


def sort_dict_by_keys(input_dict):
    sorted_keys = sorted(input_dict.keys(), key=lambda k: k.split("->")[-1].split("(")[0])
    sorted_dict = {k: input_dict[k] for k in sorted_keys}
    return dict(
        sorted(
            sorted_dict.items(),
            key=lambda x: [int(n) for n in re.findall(r"\((\d+)\)", x[0])],
            reverse=True,
        )
    )


def summarize_event_level_loss(records, node_count_mode="parsimony"):
    """Build event-level loss-node and loss-pattern labels for export."""
    event_level = {}
    for rec in records:
        event_key = (
            str(rec.get("tre_id", LOSS_TYPE_NA)),
            str(rec.get("gd_id", LOSS_TYPE_NA)),
            str(rec.get("gd_level_name", LOSS_TYPE_NA)),
        )
        entry = event_level.setdefault(
            event_key,
            {
                "nodes": [],
                "node_seen": set(),
                "pairs": [],
                "pair_seen": set(),
            },
        )

        loss_node = infer_loss_node_from_path(rec.get("loss_path", LOSS_TYPE_NA))
        if (
            loss_node
            and loss_node != LOSS_TYPE_NA
            and re.fullmatch(r"S\d+", loss_node)
            and loss_node not in entry["node_seen"]
        ):
            entry["node_seen"].add(loss_node)
            entry["nodes"].append(loss_node)

        for node_name, loss_type in rec.get("effective_transition_keys_for_export", []) or []:
            if not re.fullmatch(r"S\d+", str(node_name)):
                continue
            pair = (str(node_name), str(loss_type))
            if pair in entry["pair_seen"]:
                continue
            entry["pair_seen"].add(pair)
            entry["pairs"].append(pair)

    for entry in event_level.values():
        entry["loss_node_text"] = ";".join(entry["nodes"]) if entry["nodes"] else LOSS_TYPE_NA
        entry["loss_pattern_text"] = (
            ";".join(f"{node_name}:{loss_type}" for node_name, loss_type in entry["pairs"])
            if entry["pairs"]
            else LOSS_TYPE_NA
        )
    return event_level


def write_gd_loss_csv(records, output_file="gd_loss.csv", numed_sptree_path="numed_sptree.nwk", node_count_mode="parsimony"):
    """Write the single standard GD-loss output table.

    The exported CSV also keeps the legacy path-expanded node/species columns so
    users can trace one loss path directly from the same table.
    """
    if not records:
        logger.warning("No GD-loss records found. Creating placeholder CSV.")
        pd.DataFrame(
            {"Message": ["No gene loss paths were detected or matched the specified filters."]}
        ).to_csv(output_file, index=False)
        logger.info("Placeholder CSV saved as %s", output_file)
        return

    internal_node_columns = []
    species_columns = []
    if os.path.exists(numed_sptree_path):
        try:
            numbered_tree = PhyloTree(numed_sptree_path, format=1)
            internal_node_columns = [
                node.name for node in numbered_tree.traverse("preorder")
                if not node.is_leaf() and node.name
            ]
            species_columns = [leaf.name for leaf in numbered_tree.iter_leaves()]
            internal_node_columns = sorted(
                set(internal_node_columns),
                key=node_sort_key,
                reverse=True,
            )
        except Exception as exc:
            logger.warning("Failed to parse %s for path-expanded columns: %s", numed_sptree_path, exc)

    if not internal_node_columns or not species_columns:
        all_nodes = []
        all_species = []
        for rec in records:
            parsed = parse_loss_path(rec.get("loss_path", LOSS_TYPE_NA))
            for name, _ in parsed:
                if re.fullmatch(r"S\d+", name):
                    all_nodes.append(name)
                else:
                    all_species.append(name)
        internal_node_columns = sorted(set(all_nodes), key=node_sort_key, reverse=True)
        species_columns = sorted(set(all_species))

    all_path_columns = internal_node_columns + species_columns

    event_level_loss = summarize_event_level_loss(records, node_count_mode=node_count_mode)

    output_rows = []
    for rec in records:
        loss_path = rec.get("loss_path", LOSS_TYPE_NA)
        parsed = parse_loss_path(loss_path)
        step_map = {name: str(count) for name, count in parsed}

        tre_id = str(rec.get("tre_id", LOSS_TYPE_NA))
        gd_id = str(rec.get("gd_id", LOSS_TYPE_NA))
        gd_level_name = str(rec.get("gd_level_name", LOSS_TYPE_NA))
        event_key = (tre_id, gd_id, gd_level_name)
        event_info = event_level_loss.get(
            event_key,
            {"loss_node_text": LOSS_TYPE_NA, "loss_pattern_text": LOSS_TYPE_NA},
        )

        output_row = {
            "Tree ID": tre_id,
            "GD ID": gd_id,
            "GD burst node": gd_level_name,
            "GD loss node": event_info["loss_node_text"],
            "GD loss pattern": event_info["loss_pattern_text"],
            "GD support": rec.get("support", LOSS_TYPE_NA),
            "Species": rec.get("species_display", LOSS_TYPE_NA),
            "GD subclade A genes": rec.get("gene1", LOSS_TYPE_NA),
            "GD subclade B genes": rec.get("gene2", LOSS_TYPE_NA),
            "GD subclade A gene count": rec.get("left_gene_n", 0),
            "GD subclade B gene count": rec.get("right_gene_n", 0),
            "Species final loss type": rec.get("loss_type", LOSS_TYPE_NA),
        }
        for col in all_path_columns:
            output_row[col] = step_map.get(col, LOSS_TYPE_NA)
        output_row["GD loss path"] = loss_path
        output_rows.append(output_row)

    output_df = pd.DataFrame(output_rows)

    output_df["_tree_sort"] = output_df["Tree ID"].astype(str)
    output_df["_gd_sort"] = pd.to_numeric(output_df["GD ID"], errors="coerce").fillna(SORT_SENTINEL)
    output_df["_burst_sort"] = output_df["GD burst node"].map(node_sort_key)
    output_df["_species_sort"] = output_df["Species"].astype(str)
    output_df = output_df.sort_values(
        by=["_tree_sort", "_gd_sort", "_burst_sort", "_species_sort"],
        kind="stable",
    ).reset_index(drop=True)
    output_df = output_df.drop(columns=["_tree_sort", "_gd_sort", "_burst_sort", "_species_sort"])

    output_df.to_csv(output_file, index=False)
    logger.info("CSV report successfully generated: %s", output_file)


def write_gd_loss_node_summary_tsv(records, output_file="gd_loss_node_summary.tsv", node_count_mode="parsimony"):
    """Write node-level GD-loss summary for the selected counting mode.

    The table mirrors the statistics used by ``GD_Loss_Visualizer``:

    - ``GD event count`` is the number of unique GD events born at the node.
    - ``Node loss count`` and mode-specific ``P2_0/P2_1/P1_0`` or
      ``C2_0/C2_1/C1_0`` summarize node-level loss transitions. In accumulate
      mode, every species path transition is counted. In parsimony mode, shared
      events are collapsed by
      ``Tree ID + GD ID + Node + Type``.
    """
    mode = str(node_count_mode or "parsimony").strip().lower()
    if mode not in {"parsimony", "accumulate"}:
        mode = "parsimony"

    if not records:
        prefix = "P" if str(node_count_mode or "parsimony").strip().lower() == "parsimony" else "C"
        pd.DataFrame(
            columns=[
                "Node", "Mode", "GD event count", "Node loss count",
                f"{prefix}2_0", f"{prefix}2_1", f"{prefix}1_0",
            ]
        ).to_csv(output_file, sep="\t", index=False)
        logger.info("Node summary TSV generated: %s", output_file)
        return

    gd_birth_sets = defaultdict(set)
    node_loss_counts = defaultdict(lambda: {LOSS_TYPE_2_0: 0, LOSS_TYPE_2_1: 0, LOSS_TYPE_1_0: 0})
    seen_parsimony_events = set()

    for rec in records:
        tree_id = str(rec.get("tre_id", LOSS_TYPE_NA))
        gd_id = str(rec.get("gd_id", LOSS_TYPE_NA))
        event_key = (tree_id, gd_id)
        birth_node = str(rec.get("gd_level_name", LOSS_TYPE_NA))
        if birth_node and birth_node != LOSS_TYPE_NA:
            gd_birth_sets[birth_node].add(event_key)

        events = rec.get("effective_transition_keys_for_export", []) or []
        for node_name, transition_type in events:
            transition_type = str(transition_type).replace("-", "_")
            if transition_type not in {LOSS_TYPE_2_0, LOSS_TYPE_2_1, LOSS_TYPE_1_0}:
                continue
            if mode == "parsimony":
                logical_event = (tree_id, gd_id, node_name, transition_type)
                if logical_event in seen_parsimony_events:
                    continue
                seen_parsimony_events.add(logical_event)
            node_loss_counts[node_name][transition_type] += 1

    all_nodes = sorted(
        set(gd_birth_sets.keys()) | set(node_loss_counts.keys()),
        key=node_sort_key,
        reverse=True,
    )

    rows = []
    count_prefix = "P" if mode == "parsimony" else "C"
    for node in all_nodes:
        c20 = node_loss_counts[node][LOSS_TYPE_2_0]
        c21 = node_loss_counts[node][LOSS_TYPE_2_1]
        c10 = node_loss_counts[node][LOSS_TYPE_1_0]
        rows.append(
            {
                "Node": node,
                "Mode": mode,
                "GD event count": len(gd_birth_sets.get(node, set())),
                "Node loss count": c20 + c21 + c10,
                f"{count_prefix}2_0": c20,
                f"{count_prefix}2_1": c21,
                f"{count_prefix}1_0": c10,
            }
        )

    pd.DataFrame(rows).to_csv(output_file, sep="\t", index=False)
    logger.info("Node summary TSV generated: %s", output_file)



# ======================================================
# Section 7: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="GD Loss Tracker")
    parser.add_argument("--sptree", required=True, help="Species tree path")
    parser.add_argument("--gf", required=True, help="Gene family list path")
    args = parser.parse_args()
    out = "outfile"
    sptree = PhyloTree(args.sptree)
    num_sptree(sptree)
    tre_dic = read_and_return_dict(args.gf)

    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, taxa2voucher_dic = gene_id_transfer(
        tre_dic
    )

    os.makedirs(out, exist_ok=True)
    sp_dic, path2_treeid_dic, filtered_records = get_path_str_num_dic(
        tre_dic, sptree,
        gene2new_named_gene_dic, new_named_gene2gene_dic,
        voucher2taxa_dic, taxa2voucher_dic,
    )
    write_gd_loss_csv(filtered_records, "gd_loss.csv", "numed_sptree.nwk", node_count_mode="parsimony")
    write_gd_loss_node_summary_tsv(filtered_records, "gd_loss_node_summary.tsv", node_count_mode="parsimony")

    split_dicts = split_dict_by_first_last_char(sp_dic)
