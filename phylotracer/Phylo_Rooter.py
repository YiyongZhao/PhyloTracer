
"""
Rooting utilities and RF-guided selection for the PhyloTracer pipeline.

This module generates candidate rooted trees, evaluates them using depth and
duplication statistics, and selects optimal roots with RF-based refinement.
"""

import logging
import math
import os
import shutil

logger = logging.getLogger(__name__)
from itertools import combinations

import numpy as np
import pandas as pd
from ete3 import PhyloTree
from tqdm import tqdm

from phylotracer import (
    annotate_gene_tree,
    calculate_species_num,
    compute_tip_to_root_branch_length_variance,
    find_dup_node,
    gene_id_transfer,
    get_species_list,
    get_species_set,
    is_rooted,
    read_and_return_dict,
    read_phylo_tree,
    rename_input_tre,
    root_tre_with_midpoint_outgroup,
)
from phylotracer.BranchLength_NumericConverter import write_tree_to_newick
from phylotracer.Ortho_Retriever import (
    extract_tree,
    iterator,
    offcut_tre,
    rename_OGs_tre_name,
)

# --------------------------
# 1. Helper Functions
# --------------------------


def rename_output_tre(tree: object, name_mapping: dict, tree_id: str, output_dir: str) -> None:
    """Rename tree nodes and save to Newick format.

    Args:
        tree (object): Phylogenetic tree object.
        name_mapping (dict): Mapping from old names to new names.
        tree_id (str): Identifier for the tree.
        output_dir (str): Directory to save the Newick file.

    Returns:
        None

    Assumptions:
        Tree nodes have names that can be mapped by ``name_mapping``.
    """
    # Restore names
    for node in tree.traverse():
        if node.name in name_mapping:
            node.name = name_mapping[node.name]
    tree_str = tree.write(format=0)
    write_tree_to_newick(tree_str, tree_id, output_dir)
    

def get_species_map_and_depth(species_tree: object) -> dict:
    """Pre-compute species-tree node depths for fast lookup.

    Args:
        species_tree: ete3 Tree representing the species phylogeny.

    Returns:
        dict: Mapping from tree nodes to their depth from the root.

    Assumptions:
        The tree is connected and rooted.
    """
    species_depth = {}
    for node in species_tree.traverse("preorder"):
        if node.is_root():
            depth = 0
        else:
            depth = species_depth[node.up] + 1
        species_depth[node] = depth
    return species_depth


def annotate_mapped_depths(gene_tree: object, species_tree: object) -> object:
    """Annotate gene-tree nodes with mapped species-tree depths.

    Args:
        gene_tree (object): Gene tree with species-labeled leaves.
        species_tree (object): Species tree for LCA mapping.

    Returns:
        object: The same gene tree with ``mapped_depth`` features attached.

    Assumptions:
        Leaf names encode species identifiers compatible with the species tree.
    """
    if not gene_tree or not species_tree:
        return gene_tree

    sp_depths = get_species_map_and_depth(species_tree)

    cache = {}
    for node in gene_tree.traverse():
        species_set = get_species_set(node)
        if not species_set:
            node.add_feature("mapped_depth", 0)
            continue

        key = frozenset(species_set)
        if key in cache:
            node.add_feature("mapped_depth", cache[key])
            continue

        try:
            if len(species_set) == 1:
                species_name = list(species_set)[0]
                mapped_node = species_tree & species_name
            else:
                mapped_node = species_tree.get_common_ancestor(species_set)

            depth = sp_depths.get(mapped_node, 0)

            cache[key] = depth
            node.add_feature("mapped_depth", depth)
        except Exception:
            node.add_feature("mapped_depth", 0)
    return gene_tree


def get_species_tree_basal_set(species_tree_obj: object) -> set:
    """Select a basal species set from the species tree.

    Args:
        species_tree_obj (object): Rooted species tree.

    Returns:
        set: Species names in the smaller child clade of the root.

    Assumptions:
        The species tree root has at most two children.
    """
    if len(species_tree_obj.get_children()) < 2:
        return set(species_tree_obj.get_leaf_names())
    children = species_tree_obj.get_children()
    child_0_leaves = set(children[0].get_leaf_names())
    child_1_leaves = set(children[1].get_leaf_names())
    return child_0_leaves if len(child_0_leaves) < len(child_1_leaves) else child_1_leaves


def get_dynamic_basal_set(gene_tree_species_set: set, species_tree_root: object) -> set:
    """Compute a dynamic basal set based on gene-tree species coverage.

    Args:
        gene_tree_species_set (set): Species present in the gene tree.
        species_tree_root (object): Root node of the species tree.

    Returns:
        set: Basal species candidates for the current gene tree.

    Assumptions:
        The species tree is bifurcating; behavior for polytomies follows the
        original recursive logic.
    """
    # Recursion termination: reached a leaf node
    if species_tree_root.is_leaf():
        return set([species_tree_root.name])

    children = species_tree_root.get_children()

    # Assuming the species tree is bifurcating (if trifurcating, logic needs adjustment).
    # If there are multiple children, check which divergence the gene tree spans.

    # Get all species names under the two child branches of the current species tree node.
    # Note: Caching could optimize this, but species trees are usually small.
    clade_sets = []
    for child in children:
        clade_sets.append(set(child.get_leaf_names()))

    # Check which clades contain species from the gene tree
    present_flags = []
    for c_set in clade_sets:
        if not gene_tree_species_set.isdisjoint(c_set):
            present_flags.append(True)
        else:
            present_flags.append(False)

    # --- Core Logic ---

    # Case A: Gene tree species span the current divergence point (present in at least two sub-clades)
    if sum(present_flags) >= 2:
        # This is the earliest divergence point we are looking for.
        # We select the clade with fewer leaves as the "relative outgroup".

        valid_children = []
        for i, has_gene in enumerate(present_flags):
            if has_gene:
                valid_children.append((children[i], len(clade_sets[i])))

        # Sort by clade size in the species tree (ascending)
        valid_children.sort(key=lambda x: x[1])

        # Return all species in the smallest clade as the Basal Set
        best_outgroup_node = valid_children[0][0]
        return set(best_outgroup_node.get_leaf_names())

    # Case B: Gene tree species fall entirely within one sub-clade
    elif sum(present_flags) == 1:
        for i, has_gene in enumerate(present_flags):
            if has_gene:
                return get_dynamic_basal_set(gene_tree_species_set, children[i])

    # Case C: No gene tree species found (should not happen theoretically), return empty
    return set()


def get_all_rerooted_trees_filtered(tree: object, basal_species_set: set, newid2oldid: dict) -> list:
    """Generate candidate rerooted trees using basal taxa and MRCA strategies.

    Args:
        tree (object): Gene tree in voucher-based naming.
        basal_species_set (set): Species identifiers defining basal candidates.
        newid2oldid (dict): Mapping from renamed genes to original names for logging.

    Returns:
        list: List of candidate rerooted trees.

    Assumptions:
        Leaf names encode species as the prefix before ``_``.
    """
    rerooted_trees = []

    # To deduplicate added MRCA nodes (prevent duplicate topologies)
    # Note: We store the set of leaf names composing the clade to check for duplicate topologies
    # because node objects in copies are different.
    added_topologies = set()

    # === Collect all leaves belonging to outgroup genes (from original tree) ===
    all_basal_leaves = []
    for leaf in tree.get_leaves():
        leaf_species = leaf.name.split('_')[0] if '_' in leaf.name else leaf.name
        if leaf_species in basal_species_set:
            all_basal_leaves.append(leaf)

    if not all_basal_leaves:
        logger.info("No basal leaves found. Skipping rooting.")
        return rerooted_trees

    # old_names = [newid2oldid.get(l.name, l.name) for l in all_basal_leaves]
    # print(f"[Rooting] Found {len(all_basal_leaves)} basal leaves: {old_names}")
    is_monophyletic = False
    try:
        if len(all_basal_leaves) == 1:
            # Only one gene, inherently monophyletic
            is_monophyletic = True
        else:
            # Get the Most Recent Common Ancestor (MRCA) of all outgroups on the current tree
            mrca_node = tree.get_common_ancestor(all_basal_leaves)

            # Get all actual leaf names under this MRCA node
            mrca_leaves_names = set(mrca_node.get_leaf_names())
            # Get the target outgroup leaf names we selected
            basal_leaves_names = set([leaf.name for leaf in all_basal_leaves])

            # Criterion: Leaves under MRCA must exactly match selected outgroup leaves
            # If inner group genes are present in mrca_leaves_names, sets won't match, indicating non-monophyly
            if mrca_leaves_names == basal_leaves_names:
                is_monophyletic = True
    except Exception as e:
        logger.warning("Monophyly check failed: %s", e)
        is_monophyletic = False

    if is_monophyletic:
        # --- Strategy A: Optimization Mode (Single species/multi-copy clustered -> Root by MRCA directly) ---
        # print("[Rooting] Monophyletic outgroup detected. Optimizing: Rooting by MRCA only.")
        try:
            tree_copy = tree.copy()

            # Find equivalent node in the new tree
            basal_names = [leaf.name for leaf in all_basal_leaves]

            if len(basal_names) == 1:
                # Case with only one gene
                target_node = tree_copy & basal_names[0]
                tree_copy.set_outgroup(target_node)
                rerooted_trees.append(tree_copy)
            else:
                # Case with multiple copies but monophyletic
                nodes_in_copy = [tree_copy & name for name in basal_names]
                full_mrca_copy = tree_copy.get_common_ancestor(nodes_in_copy)

                # Root only if MRCA is not already the tree root
                if full_mrca_copy != tree_copy:
                    tree_copy.set_outgroup(full_mrca_copy)
                    rerooted_trees.append(tree_copy)

        except Exception as e:
            logger.warning("Failed optimized MRCA rooting: %s", e)

    else:
        # --- Strategy B: Comprehensive Mode (Outgroup is scattered -> Exhaustive search) ---
        # print("[Rooting] Outgroup is scattered/paraphyletic. Using exhaustive search.")
        successful_single = 0
        successful_pairwise = 0

        # === Part 1: Root by each individual outgroup gene ===
        for leaf in all_basal_leaves:
            # 1. Copy the tree
            tree_copy = tree.copy()

            # 2. [FIX] Find the equivalent node in the COPIED tree
            # Use the 'search_nodes' or simple '&' operator if names are unique
            try:
                target_node = tree_copy & leaf.name

                # 3. Set outgroup on the copied tree using the copied node
                tree_copy.set_outgroup(target_node)
                rerooted_trees.append(tree_copy)
                successful_single += 1
            except Exception as e:
                old_name = newid2oldid.get(leaf.name, leaf.name)
                logger.warning("Failed single-leaf rooting with %s: %s", old_name, e)

        # === Part 2: Root by MRCA of pairwise combinations of outgroup genes ===
        if len(all_basal_leaves) >= 2:
            for leaf1, leaf2 in combinations(all_basal_leaves, 2):
                try:
                    tree_copy = tree.copy()

                    # [FIX] Find equivalent nodes in the COPY
                    l1_copy = tree_copy & leaf1.name
                    l2_copy = tree_copy & leaf2.name

                    # Get MRCA in the COPY
                    pair_mrca = tree_copy.get_common_ancestor([l1_copy, l2_copy])

                    # Deduplication logic (based on topology/leaf names)
                    # Using sorted leaf names under the MRCA as a signature
                    mrca_leaves_signature = frozenset(pair_mrca.get_leaf_names())

                    if mrca_leaves_signature in added_topologies:
                        continue
                    added_topologies.add(mrca_leaves_signature)

                    # [FIX] Check if MRCA is already the root
                    if pair_mrca != tree_copy:
                        tree_copy.set_outgroup(pair_mrca)

                    rerooted_trees.append(tree_copy)
                    successful_pairwise += 1

                except Exception as e:
                    old1 = newid2oldid.get(leaf1.name, leaf1.name)
                    old2 = newid2oldid.get(leaf2.name, leaf2.name)
                    logger.warning("Failed pairwise MRCA rooting for %s & %s: %s", old1, old2, e)

        # Optional: Root by the MRCA of all outgroup genes
        try:
            tree_copy = tree.copy()
            # [FIX] Collect all corresponding leaves in the COPY
            basal_leaves_names = [leaf.name for leaf in all_basal_leaves]
            # Find MRCA in the COPY
            # We need to find nodes in tree_copy that match these names
            nodes_in_copy = [tree_copy & name for name in basal_leaves_names]

            full_mrca = tree_copy.get_common_ancestor(nodes_in_copy)

            mrca_leaves_signature = frozenset(full_mrca.get_leaf_names())

            if mrca_leaves_signature not in added_topologies:
                # [FIX] Check if MRCA is already the root
                if full_mrca != tree_copy:
                    tree_copy.set_outgroup(full_mrca)

                rerooted_trees.append(tree_copy)
                added_topologies.add(mrca_leaves_signature)
        except Exception as e:
            logger.warning("Failed full basal MRCA rooting: %s", e)

        # print(f"[Rooting] Generated {successful_single} single-leaf + {successful_pairwise} pairwise MRCA candidates "
        #     f"(total {len(rerooted_trees)} unique candidates).")

    return rerooted_trees


def calculate_tree_statistics(tree: object, species_tree: object) -> tuple:
    """Compute fast, RF-free statistics for candidate trees.

    Args:
        tree (object): Candidate rooted tree.
        species_tree (object): Species tree used for GD and overlap metrics.

    Returns:
        tuple: (deep, var, GD, species_overlap) summary statistics.

    Assumptions:
        The tree is bifurcating at the root.
    """
    if len(tree.children) < 2:
        return 0, 0, 0, 0
    up_clade = tree.children[1]
    down_clade = tree.children[0]

    var_up = compute_tip_to_root_branch_length_variance(up_clade)
    var_down = compute_tip_to_root_branch_length_variance(down_clade)
    var = abs(var_up - var_down)

    # Use pre-calculated mapped_depth
    if len(up_clade.get_leaf_names()) > len(down_clade.get_leaf_names()):
        deep = getattr(down_clade, 'mapped_depth', 0)
    else:
        deep = getattr(up_clade, 'mapped_depth', 0)

    species_overlap, GD = calculate_species_overlap_gd_num(tree, species_tree)
    return deep, var, GD, species_overlap


def calculate_species_overlap_gd_num(gene_tree: object, species_tree: object) -> float:
    """Compute species overlap ratio and GD count from duplication nodes.

    Args:
        gene_tree (object): Gene tree to analyze.
        species_tree (object): Species tree for duplication detection.

    Returns:
        float: Overlap ratio between the two child clades of the largest GD node.
        int: Number of detected duplication nodes.

    Assumptions:
        Duplication detection is performed by ``find_dup_node``.
    """
    dup_nodes = find_dup_node(gene_tree, species_tree)
    if not dup_nodes:
        return 0.0, 0
    largest_tree = max(dup_nodes, key=lambda node: len(node.get_leaves()))
    up_clade = largest_tree.children[1]
    down_clade = largest_tree.children[0]
    species_list_a = get_species_list(up_clade)
    species_list_b = get_species_list(down_clade)

    union_len = len(set(species_list_a) | set(species_list_b))
    if union_len == 0:
        return 0.0, len(dup_nodes)

    overlap_ratio = len(set(species_list_a) & set(species_list_b)) / union_len
    return overlap_ratio, len(dup_nodes)

# --------------------------
# 2. Core RF Calculation Logic
# --------------------------


def calculate_rf_strategy(
    tree: object,
    renamed_species_tree: object,
    renamed_length_dict: dict,
    gene_to_new_name: dict,
    new_name_to_gene: dict,
    tree_path: str,
    tree_id: str
) -> float:
    """Select RF calculation strategy based on copy number.

    Args:
        tree (object): Gene tree in voucher-based naming.
        renamed_species_tree (object): Species tree with renamed taxa.
        renamed_length_dict (dict): Length dictionary for offcut processing.
        gene_to_new_name (dict): Mapping from original gene IDs to renamed IDs.
        new_name_to_gene (dict): Mapping from renamed IDs to original IDs.
        tree_path (str): Path to the original gene tree file.
        tree_id (str): Tree identifier.

    Returns:
        float: RF distance score accumulated across relevant subtrees.

    Assumptions:
        Single-copy trees are compared directly; multi-copy trees are split into
        principal and minor ortholog sets using existing helper functions.
    """
    species_list = get_species_list(tree)
    species_set = get_species_set(tree)

    # --- Case A: Single Copy ---
    if len(species_list) == len(species_set):
        return calculate_RF_distance(tree, renamed_species_tree)

    # --- Case B: Multi Copy (Integrated User Logic) ---
    else:
        # 1. Split Principal Gene Set and Offcuts
        principal_gene_set, filtered_offcut_ev_seqs = offcut_tre(tree, renamed_length_dict)
        minor_orthologs = []
        minor_orthologs = iterator(filtered_offcut_ev_seqs, tree, gene_to_new_name, minor_orthologs, tree_path, renamed_length_dict)
        ordered_name_OG_list = rename_OGs_tre_name(principal_gene_set, minor_orthologs, tree_id)
        RF = 0
        for tre_name, OG_set in ordered_name_OG_list:
            OG_list = [new_name_to_gene[OG] for OG in OG_set]
            phylo_tree_0 = read_phylo_tree(tree_path)
            phylo_tree = root_tre_with_midpoint_outgroup(phylo_tree_0)
            phylo_tree_OG_list = extract_tree(OG_list, phylo_tree)
            phylo_tree_OG_list = rename_input_tre(phylo_tree_OG_list, gene_to_new_name)
            RF += calculate_RF_distance(phylo_tree_OG_list, renamed_species_tree)

        return RF


def calculate_RF_distance(Phylo_t_OG_L: object, sptree: object) -> int:
    """Compute Robinson-Foulds distance with basic preprocessing.

    Args:
        Phylo_t_OG_L (object): Gene tree to compare.
        sptree (object): Species tree reference.

    Returns:
        int: Robinson-Foulds distance; defaults to 100 on error.

    Assumptions:
        Leaf names are formatted as ``Species_Gene`` and are reduced to species IDs.
    """
    tcopy = Phylo_t_OG_L.copy()
    for leaf in tcopy:
        if "_" in leaf.name:
            leaf.name = leaf.name.split('_')[0]
    try:
        rf, _ = tcopy.robinson_foulds(sptree)[:2]
        return rf
    except Exception:
        return 100  # Penalty for error

# --------------------------
# 3. Scoring and Normalization
# --------------------------


def normalize_and_score(df: pd.DataFrame, weights: dict, include_rf=True):
    """Normalize metrics and compute a composite score.

    Args:
        df: pandas DataFrame with metric columns.
        weights (dict): Weighting scheme for metrics.
        include_rf (bool): Whether to include RF in the score.

    Returns:
        pandas Series: Score values for each row.

    Assumptions:
        Metrics are numeric and higher overlap is better while deep/var/GD are lower.
    """
    def norm(values: pd.Series) -> pd.Series:
        """Normalize a series of values to [0, 1] range.

        Args:
            values (pd.Series): Input metric values.

        Returns:
            pd.Series: Normalized values.
        """
        if values.max() != values.min():
            return (values - values.min()) / (values.max() - values.min())
        return pd.Series(0.0, index=values.index)

    # deep/var/GD: smaller is better; overlap: larger is better
    norm_deep = norm(df["deep"])
    norm_var = norm(df["var"])
    norm_GD = norm(df["GD"])
    norm_species_overlap = norm(df["species_overlap"])

    weighted_norm_RF = 0
    if include_rf and "RF" in df.columns:
        norm_RF = norm(df["RF"])
        weighted_norm_RF = norm_RF * weights.get("RF", 0)

    weighted_norm_overlap = norm_species_overlap * weights.get("species_overlap", 0)
    score = (
        norm_deep * weights.get("deep", 0) +
        norm_var * weights.get("var", 0) +
        norm_GD * weights.get("GD", 0) -
        weighted_norm_overlap +
        weighted_norm_RF
    )
    return score

# --------------------------
# 4. Main Process
# --------------------------


def root_main(
    tree_dict: dict,
    gene_to_new_name: dict,
    renamed_length_dict: dict,
    new_name_to_gene: dict,
    renamed_species_tree: object
) -> None:
    """Root gene trees and select optimal candidates using staged scoring.

    Args:
        tree_dict (dict): Mapping from tree identifiers to tree file paths.
        gene_to_new_name (dict): Mapping from original gene IDs to renamed IDs.
        renamed_length_dict (dict): Length dictionary used in offcut processing.
        new_name_to_gene (dict): Mapping from renamed IDs to original IDs.
        renamed_species_tree (object): Species tree with renamed taxa.

    Returns:
        None

    Assumptions:
        Input trees are readable Newick files and species labels are consistent.
    """
    cwd = os.getcwd()
    default_dir = "rooted_trees"
    dir_path = (
        cwd
        if os.path.basename(os.path.normpath(cwd)) == default_dir
        else os.path.join(cwd, f"{default_dir}/")
    )
    if dir_path != cwd and os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path, exist_ok=True)

    # Pre-processing: Get basal taxa
    basal_species_set = get_species_tree_basal_set(renamed_species_tree)
    logger.info("Basal Filter: %d species", len(basal_species_set))

    pbar = tqdm(total=len(tree_dict), desc="Processing trees", unit="tree")
    stat_matrix = []

    stage1_weights_multi = {"deep": 0.3, "var": 0.1, "RF": 0.0, "GD": 0.5, "species_overlap": 0.1}
    stage1_weights_single = {"deep": 0.7, "var": 0.3, "RF": 0.0}

    try:
        for tree_id, tree_path in tree_dict.items():
            pbar.set_description(f"Processing {tree_id}")

            # 1. Read and initial processing
            Phylo_t0 = read_phylo_tree(tree_path)
            Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
            Phylo_t2 = rename_input_tre(Phylo_t1, gene_to_new_name)

            # 2. Handle simple cases
            if len(get_species_set(Phylo_t2)) == 1 or len(get_species_list(Phylo_t2)) <= 3:
                rename_output_tre(Phylo_t2, new_name_to_gene, tree_id, dir_path)
                pbar.update(1)
                continue

            # 3. Depth annotation & Basal filtering
            Phylo_t2 = annotate_mapped_depths(Phylo_t2, renamed_species_tree)

            current_gene_species = set(get_species_list(Phylo_t2))
            dynamic_basal_set = get_dynamic_basal_set(current_gene_species, renamed_species_tree)
            root_list = get_all_rerooted_trees_filtered(Phylo_t2, dynamic_basal_set, new_name_to_gene)

            if not root_list:
                logger.warning("No valid root candidates after basal filter.")
                rename_output_tre(Phylo_t2, new_name_to_gene, tree_id, dir_path)
                pbar.update(1)
                continue

            # ==========================================
            # Stage 1: Fast Screening (No RF)
            # ==========================================
            tree_objects = {}
            temp_stats = []

            for n, tree in enumerate(root_list):
                # tree is already in renamed Taxa ID format
                tree.resolve_polytomy(recursive=True)
                tree_key = f"{tree_id}_{n+1}"
                tree_objects[tree_key] = tree
                # print(rename_input_tre(tree, new_name_to_gene))

                deep, var, GD, species_overlap = calculate_tree_statistics(tree, renamed_species_tree)
                # print(f"deep: {deep}, var: {var}, GD: {GD}, species_overlap: {species_overlap}")
                # print('='*50)
                temp_stats.append({
                    "Tree": tree_key,
                    "deep": deep, "var": var, "GD": GD,
                    "species_overlap": species_overlap, "RF": 0,
                    "tree_obj_ref": tree  # Store object reference for Stage 2
                })

            current_df = pd.DataFrame(temp_stats)

            # Weight selection
            is_multi_copy = len(get_species_list(Phylo_t2)) != len(get_species_set(Phylo_t2))
            used_weights = stage1_weights_multi if is_multi_copy else stage1_weights_single

            # Initial scoring
            current_df["score"] = normalize_and_score(current_df, used_weights, include_rf=False)

            # ==========================================
            # Stage 2: Fine Screening (Calculate RF for Top N)
            # ==========================================
            total_candidates = len(current_df)
            # Strategy: Top 40% or Top 20 (Keep the funnel loose to avoid missing trees with good RF)
            top_n = max(20, math.ceil(total_candidates * 0.8))
            top_n = min(top_n, total_candidates)

            # 1. Select Top N with best initial scores
            top_candidates_df = current_df.nsmallest(top_n, "score").copy()

            # 2. Calculate RF for Top N
            for idx, row in top_candidates_df.iterrows():
                target_tree = row["tree_obj_ref"]

                # Calculate RF
                rf_val = calculate_rf_strategy(
                    tree=target_tree,
                    renamed_species_tree=renamed_species_tree,
                    renamed_length_dict=renamed_length_dict,
                    gene_to_new_name=gene_to_new_name,
                    new_name_to_gene=new_name_to_gene,
                    tree_path=tree_path,
                    tree_id=tree_id
                )

                top_candidates_df.at[idx, "RF"] = rf_val

            # 3. [Core Modification] Sort directly to select the best, no re-weighting
            # Logic: Prioritize RF (smaller is better), if RF is equal, check Stage 1 Score (smaller is better)
            best_row = top_candidates_df.sort_values(
                by=["RF", "score"],
                ascending=[True, True]
            ).iloc[0]

            # ==========================================
            # Output
            # ==========================================
            best_tree_key = best_row["Tree"]
            best_tree = tree_objects[best_tree_key]

            rename_output_tre(best_tree, new_name_to_gene, tree_id, dir_path)

            record = best_row.to_dict()
            del record["tree_obj_ref"]
            stat_matrix.append(record)

            pbar.update(1)

        # Save report
        stat_df = pd.DataFrame(stat_matrix)
        if not stat_df.empty:
            stat_df["tree_id"] = stat_df["Tree"].apply(lambda x: "_".join(x.split("_")[:-1]))
            cols = ["Tree", "score", "deep", "var", "GD", "species_overlap", "RF"]
            stat_df = stat_df.sort_values(by=["tree_id", "score"])[cols]
            stat_df.to_csv("stat_matrix.csv", index=False)

    finally:
        pbar.close()


