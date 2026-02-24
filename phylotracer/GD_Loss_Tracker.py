"""
Gene duplication loss tracking for the PhyloTracer pipeline.

This module evaluates loss events following duplication nodes, summarizes
loss paths, and generates tabular reports for downstream analyses.
"""

import logging
import re
from typing import Optional

import os

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
)

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
    2-2/2-1/2-0 classification and should not be used in new code.
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
        one of ``2-2``, ``2-1``, ``2-0``, or ``NA``.
    """
    children = dup_node.get_children()
    if mapped_node is None or len(children) != 2:
        return "NA", False, False, "invalid_node"

    if species_voucher not in set(mapped_node.get_leaf_names()):
        return "NA", False, False, "out_of_mapped_node"

    left_child, right_child = children
    left_species = get_species_set(left_child)
    right_species = get_species_set(right_child)

    left_has = species_voucher in left_species
    right_has = species_voucher in right_species

    if family_species_set is not None and species_voucher not in family_species_set:
        if not include_unobserved_species:
            return "NA", left_has, right_has, "unobserved_in_family"
        confidence = "unobserved_in_family_included"
    else:
        confidence = "observable"

    if left_has and right_has:
        loss_type = "2-2"
    elif left_has or right_has:
        loss_type = "2-1"
    else:
        loss_type = "2-0"

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
        ``(path_count_types, c20, c21, c10)`` where counts are path-count
        transitions (not biological copy-number states).
    """
    if not path_str or path_str == "NA":
        return "NA", 0, 0, 0

    parsed = []
    for step in path_str.split("->"):
        m = re.search(r"(.+?)\((\d+)\)$", step.strip())
        if m:
            parsed.append((m.group(1).strip(), int(m.group(2))))

    if len(parsed) < 2:
        return "NA", 0, 0, 0

    c20 = c21 = c10 = 0
    ordered_types = []
    transition_keys = []
    for i in range(1, len(parsed)):
        prev_node = parsed[i - 1][0]
        prev_copy = parsed[i - 1][1]
        curr_copy = parsed[i][1]
        if curr_copy >= prev_copy:
            continue
        if prev_copy == 2 and curr_copy == 0:
            c20 += 1
            ordered_types.append("2-0")
            transition_keys.append((prev_node, "2-0"))
        elif prev_copy == 2 and curr_copy == 1:
            c21 += 1
            ordered_types.append("2-1")
            transition_keys.append((prev_node, "2-1"))
        elif prev_copy == 1 and curr_copy == 0:
            c10 += 1
            ordered_types.append("1-0")
            transition_keys.append((prev_node, "1-0"))

    return (",".join(ordered_types) if ordered_types else "NA"), c20, c21, c10, transition_keys


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

    tips = set([j for j in sp_dic.keys() if not j.startswith("S")])
    for k, v in sp_dic.items():
        if k.startswith("S"):
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
    node_count_mode: str = "nonaccumulate",
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
    dup_nodes = find_dup_node(genetree, renamed_sptree,50,0,2,0,1)

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
        gd_counter = {"2-2": 0, "2-1": 0, "2-0": 0, "missing_data": 0}
        gd_seen_transition_keys = set()
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
            g_a = left_genes[0] if left_genes else "NA"
            g_b = right_genes[0] if right_genes else "NA"

            orig_a = new_named_gene2gene_dic.get(g_a, g_a)
            orig_b = new_named_gene2gene_dic.get(g_b, g_b)
            pretty_path = voucher_to_pretty_path.get(species_voucher, "NA")
            path_count_types, path_count_2_0, path_count_2_1, path_count_1_0, transition_keys = (
                summarize_small_loss_types_from_path(pretty_path)
            )
            transition_keys_effective = list(transition_keys)
            if node_count_mode == "nonaccumulate" and transition_keys:
                filtered_types = []
                c20 = c21 = c10 = 0
                filtered_transition_keys = []
                for node_name, t in transition_keys:
                    key = (node_name, t)
                    if key in gd_seen_transition_keys:
                        continue
                    gd_seen_transition_keys.add(key)
                    filtered_types.append(t)
                    filtered_transition_keys.append((node_name, t))
                    if t == "2-0":
                        c20 += 1
                    elif t == "2-1":
                        c21 += 1
                    elif t == "1-0":
                        c10 += 1
                path_count_types = ",".join(filtered_types) if filtered_types else "NA"
                path_count_2_0, path_count_2_1, path_count_1_0 = c20, c21, c10
                transition_keys_effective = filtered_transition_keys

            if transition_keys_effective:
                # Node-aware transition events for downstream visualization.
                path_count_node_events = ";".join(
                    f"{node_name}|{trans_type}" for node_name, trans_type in transition_keys_effective
                )
            else:
                path_count_node_events = "NA"
            if confidence == "unobserved_in_family" and not include_unobserved_species:
                major_loss_class = "missing_data"
            elif loss_type == "2-0":
                major_loss_class = "loss_two_copies"
            elif loss_type == "2-1":
                major_loss_class = "loss_one_copy"
            elif loss_type == "2-2":
                major_loss_class = "no_loss"
            else:
                major_loss_class = "NA"

            if loss_type in gd_counter:
                gd_counter[loss_type] += 1
            elif major_loss_class == "missing_data":
                gd_counter["missing_data"] += 1

            records.append(
                {
                    "tre_id": tre_id,
                    "gd_id": gd_id,
                    "support": dup_node.support,
                    "gd_level_name": gd_level_name,
                    "species_display": voucher2taxa_dic.get(
                        species_voucher,
                        species_voucher,
                    ),
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
                    "path_count_node_events": path_count_node_events,
                    "path_count_types": path_count_types,
                    "path_count_2_0": path_count_2_0,
                    "path_count_2_1": path_count_2_1,
                    "path_count_1_0": path_count_1_0,
                }
            )

        clade_species_count = len(clade_species)
        observed_species_in_clade = len(clade_species & family_species_set)
        n22 = gd_counter["2-2"]
        n21 = gd_counter["2-1"]
        n20 = gd_counter["2-0"]
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
    node_count_mode: str = "nonaccumulate",
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

    with open("gd_loss_summary.txt", "w") as out:
        out.write(
            "tree_ID\tgd_ID\tgd_support\tlevel\tspecies\tgene1\tgene2\tloss_path\t"
            "left_gene_n\tright_gene_n\tloss_type\tleft_has\tright_has\t"
            "loss_confidence\tmajor_loss_class\t"
            "path_count_node_events\tpath_count_types\tpath_count_2_0\tpath_count_2_1\tpath_count_1_0\n"
        )

        for rec in all_records:
            if target_species_list:
                if rec["species_voucher"] not in target_voucher_set:
                    continue

            if allowed_gd_species_sets:
                gd_node = renamed_sptree & rec["gd_node_name"]
                gd_leaves_voucher = gd_node.get_leaf_names()
                gd_leaves_taxa = frozenset(
                    voucher2taxa_dic.get(leaf, leaf) for leaf in gd_leaves_voucher
                )
                if gd_leaves_taxa not in allowed_gd_species_sets:
                    continue

            out.write(
                f"{rec['tre_id']}\t{rec['gd_id']}\t{rec['support']}\t"
                f"{rec['gd_level_name']}\t"
                f"{rec['species_display']}\t"
                f"{rec['gene1']}\t{rec['gene2']}\t{rec['loss_path']}\t"
                f"{rec.get('left_gene_n', 0)}\t{rec.get('right_gene_n', 0)}\t"
                f"{rec.get('loss_type', 'NA')}\t{rec.get('left_has', 0)}\t"
                f"{rec.get('right_has', 0)}\t{rec.get('loss_confidence', 'NA')}\t"
                f"{rec.get('major_loss_class', 'NA')}\t"
                f"{rec.get('path_count_node_events', 'NA')}\t"
                f"{rec.get('path_count_types', 'NA')}\t"
                f"{rec.get('path_count_2_0', 0)}\t{rec.get('path_count_2_1', 0)}\t"
                f"{rec.get('path_count_1_0', 0)}\n"
            )

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

    sorted_sp_dict = sort_dict_by_keys(path_str_num_dic)
    with open("gd_loss_count_summary.txt", "w") as f:
        f.write("GD Loss path\tGF count\n")
        processed_sp = set()
        for k, v in sorted_sp_dict.items():
            last_char = k.split("->")[-1].split("(")[0]
            converted_last_char = voucher2taxa_dic.get(last_char, last_char)
            new_k = "->".join(k.split("->")[:-1]) + "->" + k.split("->")[-1].replace(
                last_char,
                converted_last_char,
                1,
            )
            if last_char not in processed_sp:
                f.write(f"\n{new_k}\t{v}\n")
                processed_sp.add(last_char)
            else:
                f.write(f"{new_k}\t{v}\n")

    if not path_str_num_dic:
        logger.warning("No paths matched the filters. Result is empty.")

    return path_str_num_dic, path2_treeid_dic


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


def parse_text_to_excel(file_path, output_file="gd_loss.xlsx"):
    """
    Parse loss summary text and write a multi-sheet Excel report.

    Parameters
    ----------
    file_path : str
        Path to ``gd_loss_count_summary.txt``.
    output_file : str, optional
        Output Excel file path.

    Returns
    -------
    None

    Assumptions
    -----------
    The summary file follows the expected path-count format.
    """
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        logger.warning("Summary file is missing or empty. Creating placeholder Excel.")
        placeholder_df = pd.DataFrame(
            {
                "Message": [
                    "No gene loss paths were detected or matched the specified filters."
                ]
            }
        )
        placeholder_df.to_excel(output_file, sheet_name="No_Data", index=False)
        logger.info("Placeholder Excel saved as %s", output_file)
        return

    group_dict = {}
    with open(file_path, "r") as file:
        first_line = True
        for line in file:
            if first_line:
                first_line = False
                continue
            if not line.strip():
                continue
            parts = line.strip().split("->")
            if len(parts) < 2:
                continue
            start_node = parts[0].split("(")[0]
            species = parts[-1].split("(")[0]
            key = (start_node, species)
            if key not in group_dict:
                group_dict[key] = []
            group_dict[key].append(line.strip())

    if not group_dict:
        logger.warning("No valid data groups found after parsing. Creating placeholder Excel.")
        placeholder_df = pd.DataFrame(
            {"Message": ["No gene loss paths matched the filters after parsing."]}
        )
        placeholder_df.to_excel(output_file, sheet_name="No_Data", index=False)
        logger.info("Placeholder Excel saved as %s", output_file)
        return

    def parse_lines_to_df(lines):
        header_line = lines[0]
        columns = [col.split("(")[0] for col in header_line.split("->")]
        columns.append("gd_number")

        rows = []
        for line in lines:
            if not line:
                continue
            matches = re.findall(r"\((\d+)\)", line)
            count_match = re.search(r"\s+(\d+)$", line)
            if count_match and matches:
                count = count_match.group(1)
                row = matches + [count]
                if len(row) == len(columns):
                    rows.append(row)

        if not rows:
            return None

        df = pd.DataFrame(rows, columns=columns)

        descriptions = []
        for _, row in df.iterrows():
            loss_desc = []
            for i in range(len(columns) - 2):
                if int(row[columns[i]]) < 2:
                    if i > 0:
                        loss_desc.append(f"Lost after {columns[i-1]}")
                    else:
                        loss_desc.append("Lost after root")
                    break
            desc = "No duplicate lost" if not loss_desc else loss_desc[0]
            descriptions.append(desc)

        df["Rest # of duplicates"] = descriptions
        df = df[["Rest # of duplicates"] + columns]
        return df

    with pd.ExcelWriter(output_file) as writer:
        for (start, species), lines in group_dict.items():
            dic = {}
            for ln in lines:
                path, num = ln.strip().split("\t")
                dic[path] = int(num)
            sorted_dic = sort_dict_by_keys(dic)
            sorted_lines = [f"{k}\t{v}" for k, v in sorted_dic.items()]

            df = parse_lines_to_df(sorted_lines)
            if df is None or df.empty:
                continue

            sheet_name = f"{start}_{species}"[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    logger.info("Excel report successfully generated: %s", output_file)


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
    sp_dic, path2_treeid_dic = get_path_str_num_dic(
        tre_dic, sptree,
        gene2new_named_gene_dic, new_named_gene2gene_dic,
        voucher2taxa_dic, taxa2voucher_dic,
    )

    split_dicts = split_dict_by_first_last_char(sp_dic)
