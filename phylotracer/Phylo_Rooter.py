
"""
Rooting utilities and RF-guided selection for the PhyloTracer pipeline.

This module generates candidate rooted trees, evaluates them using depth and
duplication statistics, and selects optimal roots with RF-based refinement.
"""

import logging
import os
import shutil

logger = logging.getLogger(__name__)
from itertools import combinations

import numpy as np
import pandas as pd
from ete3 import PhyloTree, Tree
from tqdm import tqdm

from phylotracer import (
    annotate_gene_tree,
    compute_tip_to_root_branch_length_variance,
    find_dup_node,
    get_species_list,
    get_species_set,
    map_species_set_to_node,
    read_phylo_tree,
    rename_input_tre,
    root_tre_with_midpoint_outgroup,
    serialize_tree_by_input_branch_length_style,
)
from phylotracer.MulRF_Distance import compute_mulrf
from phylotracer.BranchLength_NumericConverter import write_tree_to_newick

# --------------------------
# 1. Helper Functions
# --------------------------


def rename_output_tre(
    tree: object,
    name_mapping: dict,
    tree_id: str,
    output_dir: str,
    source_tree_path: str = None,
) -> None:
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
    tree_str = serialize_tree_by_input_branch_length_style(
        tree,
        source_tree_path=source_tree_path,
        fmt=0,
    )
    write_tree_to_newick(tree_str, tree_id, output_dir)
    

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


def get_species_map_and_depth(species_tree: object) -> dict:
    """Compute topological depth for every node in a species tree.

    Args:
        species_tree (object): ETE species tree.

    Returns:
        dict: Mapping from each tree node to its topological depth
        (root = 0, children of root = 1, …).
    """
    depths = {}
    for node in species_tree.traverse("levelorder"):
        if node.is_root():
            depths[node] = 0
        else:
            depths[node] = depths[node.up] + 1
    return depths


def annotate_mapped_depths(gene_tree, species_tree):
    """Annotate every gene-tree node with its mapped species-tree depth.

    Each node receives a ``mapped_depth`` feature equal to the topological
    depth of the corresponding species-tree node (determined via the
    species prefix of its descendant leaves).

    Args:
        gene_tree: ETE gene tree, or ``None``.
        species_tree: ETE species tree, or ``None``.

    Returns:
        The annotated gene tree, ``None`` when *gene_tree* is ``None``,
        or the original *gene_tree* unchanged when *species_tree* is ``None``.
    """
    if gene_tree is None:
        return None
    if species_tree is None:
        return gene_tree

    depth_map = get_species_map_and_depth(species_tree)
    sp_leaves = set(species_tree.get_leaf_names())

    for node in gene_tree.traverse():
        species = get_species_set(node)
        mapped_species = species & sp_leaves
        if not mapped_species:
            node.add_feature("mapped_depth", 0)
            continue
        try:
            if len(mapped_species) == 1:
                sp_node = species_tree & list(mapped_species)[0]
            else:
                sp_node = species_tree.get_common_ancestor(list(mapped_species))
            node.add_feature("mapped_depth", depth_map.get(sp_node, 0))
        except (ValueError, KeyError):
            node.add_feature("mapped_depth", 0)

    return gene_tree


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


def _root_signature(tree: object) -> tuple:
    """Build a stable signature of the current root split for candidate deduplication."""
    if len(tree.children) < 2:
        return tuple()
    parts = []
    for child in tree.children[:2]:
        parts.append(tuple(sorted(child.get_leaf_names())))
    return tuple(sorted(parts))


def _reroot_by_node_set(base_tree: object, node_list: list) -> list:
    """Generate rerooted candidates by setting outgroup to each node in ``node_list``."""
    rerooted = []
    seen = set()
    for node in node_list:
        try:
            leaf_names = node.get_leaf_names()
            tcopy = base_tree.copy()
            if len(leaf_names) == 1:
                target = tcopy & leaf_names[0]
            else:
                target = tcopy.get_common_ancestor(leaf_names)
            if target == tcopy:
                continue
            tcopy.set_outgroup(target)
            sig = _root_signature(tcopy)
            if sig in seen:
                continue
            seen.add(sig)
            rerooted.append(tcopy)
        except Exception as exc:
            logger.debug("Rerooting failed for node (non-fatal): %s", exc)
            continue
    return rerooted


def merge_unique_root_candidates(*candidate_groups: list) -> list:
    """Merge candidate trees and keep only unique root-split signatures."""
    merged = []
    seen = set()
    for group in candidate_groups:
        for tree in group:
            sig = _root_signature(tree)
            if sig in seen:
                continue
            seen.add(sig)
            merged.append(tree)
    return merged


def _tip_distance_mad(tree: object) -> float:
    """Compute mean absolute deviation of root-to-tip distances."""
    distances = [tree.get_distance(leaf) for leaf in tree.iter_leaves()]
    if len(distances) < 2:
        return 0.0
    median_val = float(np.median(distances))
    return float(np.mean([abs(x - median_val) for x in distances]))


def _select_statistical_root_candidates(base_tree: object, max_candidates: int = 20) -> tuple:
    """Generate MAD/MinVar root candidates from internal-node rerooting.

    Args:
        base_tree (object): Input rooted tree in renamed ID space.
        max_candidates (int): Max candidates kept per strategy.

    Returns:
        tuple: (mad_candidates, minvar_candidates)
    """
    internal_nodes = [n for n in base_tree.traverse() if (not n.is_leaf() and not n.is_root())]
    if not internal_nodes:
        return [], []

    all_candidates = _reroot_by_node_set(base_tree, internal_nodes)
    if not all_candidates:
        return [], []

    if max_candidates <= 0:
        max_candidates = 1
    keep_n = min(max_candidates, len(all_candidates))

    mad_ranked = sorted(all_candidates, key=_tip_distance_mad)
    var_ranked = sorted(all_candidates, key=compute_tip_to_root_branch_length_variance)
    return mad_ranked[:keep_n], var_ranked[:keep_n]


def calculate_tree_statistics(
    tree: object,
    species_tree: object,
) -> tuple:
    """Compute candidate-tree statistics including RF.

    Args:
        tree (object): Candidate rooted tree.
        species_tree (object): Species tree used for GD and overlap metrics.

    Returns:
        tuple: (deep, var, GD, species_overlap, gd_consistency, RF) summary statistics.

    Assumptions:
        The tree is bifurcating at the root.
    """
    if len(tree.children) < 2:
        return 0, 0, 0, 0, 0, 1.0  # RF fallback was 100, but normal range is [0,1]
    up_clade = tree.children[1]
    down_clade = tree.children[0]

    var_up = compute_tip_to_root_branch_length_variance(up_clade)
    var_down = compute_tip_to_root_branch_length_variance(down_clade)
    var = abs(var_up - var_down)

    def _mapped_topo_depth(clade: object) -> float:
        # Re-map by species set for the CURRENT rooted candidate tree.
        try:
            clade_species = get_species_set(clade)
            mapped_node = map_species_set_to_node(species_tree, clade_species)
            if mapped_node is not None:
                return float(species_tree.get_distance(mapped_node, topology_only=True) + 1)
        except Exception as exc:
            logger.debug("Species-set mapping failed (non-fatal): %s", exc)
        try:
            mapped_name = getattr(clade, "map", None)
            if mapped_name:
                try:
                    mapped_node = species_tree & mapped_name
                    return float(species_tree.get_distance(mapped_node, topology_only=True) + 1)
                except Exception as exc:
                    logger.debug("Named-map lookup failed (non-fatal): %s", exc)
        except Exception as exc:
            logger.debug("Map attribute access failed (non-fatal): %s", exc)
        fallback_depth = getattr(clade, "depth", 0)
        return float(0 if fallback_depth is None else fallback_depth)

    # Deep is defined on the side with fewer tips (proxy for outgroup side),
    # not by taking the minimum depth value across both sides.
    if len(up_clade.get_leaf_names()) > len(down_clade.get_leaf_names()):
        deep = _mapped_topo_depth(down_clade)
    else:
        deep = _mapped_topo_depth(up_clade)

    species_overlap, GD, gd_consistency = calculate_species_overlap_gd_num(tree, species_tree)

    rf_val = calculate_mulrf_rf_distance(tree, species_tree)
    return deep, var, GD, species_overlap, gd_consistency, rf_val


def calculate_species_overlap_gd_num(gene_tree: object, species_tree: object) -> tuple:
    """Compute species overlap ratio and GD count from duplication nodes.

    Args:
        gene_tree (object): Gene tree to analyze.
        species_tree (object): Species tree for duplication detection.

    Returns:
        float: Overlap ratio between the two child clades of the largest GD node.
        int: Number of detected duplication nodes.
        float: Mean GD consistency across duplication nodes.

    Assumptions:
        Duplication detection is performed by ``find_dup_node``.
    """
    dup_nodes = find_dup_node(gene_tree, species_tree)
    if not dup_nodes:
        return 0.0, 0, 0.0
    largest_tree = max(dup_nodes, key=lambda node: len(node.get_leaves()))
    if len(largest_tree.children) < 2:
        return 0.0, len(dup_nodes), 0.0
    up_clade = largest_tree.children[1]
    down_clade = largest_tree.children[0]
    species_list_a = get_species_list(up_clade)
    species_list_b = get_species_list(down_clade)

    union_len = len(set(species_list_a) | set(species_list_b))
    if union_len == 0:
        return 0.0, len(dup_nodes), 0.0

    overlap_ratio = len(set(species_list_a) & set(species_list_b)) / union_len
    gd_consistency_vals = []
    for dup_node in dup_nodes:
        if len(dup_node.children) < 2:
            continue
        left_species = set(get_species_list(dup_node.children[0]))
        right_species = set(get_species_list(dup_node.children[1]))
        if not left_species or not right_species:
            continue
        size_symmetry = min(len(left_species), len(right_species)) / max(
            len(left_species), len(right_species)
        )
        union_size = len(left_species | right_species)
        jaccard = (
            len(left_species & right_species) / union_size if union_size else 0.0
        )
        gd_consistency_vals.append(size_symmetry * jaccard)

    gd_consistency = float(np.mean(gd_consistency_vals)) if gd_consistency_vals else 0.0
    return overlap_ratio, len(dup_nodes), gd_consistency

# --------------------------
# 2. Core RF Calculation Logic
# --------------------------


def calculate_mulrf_rf_distance(gene_tree: PhyloTree, species_tree: PhyloTree, unrooted: bool = True) -> float:
    """Compute MulRF RF rate by reusing ``phylotracer.MulRF_Distance.compute_mulrf``."""
    try:
        _ = unrooted  # kept for signature compatibility
        res = compute_mulrf(gene_tree, species_tree)
        rf_rate = res.get("normalized_mulrf")
        if rf_rate is None:
            return 1.0
        return float(rf_rate)
    except Exception as exc:
        logger.debug("MulRF failed: %s", exc)
        return 1.0


# --------------------------
# 3. Scoring and Normalization
# --------------------------


def _minmax_norm(s: pd.Series, cost: bool) -> pd.Series:
    """Direction-aware min-max normalization to [0, 1], higher = better.

    Args:
        s (pd.Series): Raw metric values.
        cost (bool): True if lower raw value is better (cost metric),
            False if higher raw value is better (benefit metric).

    Returns:
        pd.Series: Values in [0, 1] where 1 means best.
            When max == min (no discriminative power), returns 0.5 (neutral).
    """
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or pd.isna(hi):  # guard against NaN in normalization
        return s
    if hi == lo:
        return pd.Series(0.5, index=s.index)
    normed = (s - lo) / (hi - lo)
    return (1.0 - normed) if cost else normed


def normalize_and_score(
    df: pd.DataFrame,
    weights: dict,
    include_rf: bool = True,
) -> pd.Series:
    """Direction-aware min-max normalization and weighted composite scoring.

    All six metrics (OD, BLV, GD, SO, GD_consistency, MulRF) are mapped to
    [0, 1] with **higher = better** before weighting, ensuring that assigned
    weights have their intended relative influence regardless of raw scale.

    Cost metrics (lower raw value is better):
        deep (OD), var (BLV), GD, RF -- normalized as (max - x)/(max - min).

    Benefit metrics (higher raw value is better):
        species_overlap (SO), gd_consistency -- normalized as (x - min)/(max - min).

    When max == min for any metric the column is set to 0.5 (neutral).

    Args:
        df (pd.DataFrame): Candidate-tree statistics with columns
            ``deep``, ``var``, ``GD``, ``species_overlap``, ``gd_consistency``,
            and optionally ``RF``.
        weights (dict): Prior-knowledge weights keyed by column name.
        include_rf (bool): Whether to include the RF column in scoring.

    Returns:
        pd.Series: Combined score per candidate. **Higher is better.**
    """
    COST_COLS    = ["deep", "var", "GD", "RF"]
    BENEFIT_COLS = ["species_overlap", "gd_consistency"]

    active_cols = ["deep", "var", "GD", "species_overlap", "gd_consistency"]
    if include_rf and "RF" in df.columns:
        active_cols.append("RF")

    # Build normalized DataFrame (all higher = better)
    normed = pd.DataFrame(index=df.index)
    for col in active_cols:
        if col not in df.columns:
            continue
        normed[col] = _minmax_norm(df[col], cost=(col in COST_COLS))

    # Re-normalize empirical weights to only the active (present) columns so
    # that dropping RF (include_rf=False) does not deflate the total score.
    raw_w = {c: weights.get(c, 0.0) for c in normed.columns}
    w_sum = sum(raw_w.values())
    if w_sum > 0:
        w = {c: v / w_sum for c, v in raw_w.items()}
    else:
        n = len(normed.columns)
        if n == 0:  # guard against empty columns
            return df
        w = {c: 1.0 / n for c in normed.columns}

    score = pd.Series(0.0, index=df.index)
    for col in normed.columns:
        score += normed[col] * w.get(col, 0.0)
    return score

# --------------------------
# 4. Main Process
# --------------------------


def root_main(
    tree_dict: dict,
    gene_to_new_name: dict,
    new_name_to_gene: dict,
    renamed_species_tree: object,
    stage1_weights: dict = None,
) -> None:
    """Root gene trees and select optimal candidates using staged scoring.

    Args:
        tree_dict (dict): Mapping from tree identifiers to tree file paths.
        gene_to_new_name (dict): Mapping from original gene IDs to renamed IDs.
        new_name_to_gene (dict): Mapping from renamed IDs to original IDs.
        renamed_species_tree (object): Species tree with renamed taxa.
        stage1_weights (dict, optional): Prior-knowledge weights for
            deep/var/GD/species_overlap/gd_consistency/RF.

    Returns:
        None

    Assumptions:
        Input trees are readable Newick files and species labels are consistent.
    """
    cwd = os.getcwd()
    explicit_output_dir = os.environ.get("PHYLTR_OUTPUT_DIR_EXPLICIT", "0") == "1"
    default_dir = "rooted_trees"
    dir_path = cwd if explicit_output_dir else (
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

    pbar = tqdm(total=len(tree_dict), desc="Processing trees", unit="tree", dynamic_ncols=True)
    stat_matrix = []

    if stage1_weights is None:
        stage1_weights = {
            "deep": 0.30,
            "var": 0.10,
            "GD": 0.30,
            "species_overlap": 0.10,
            "gd_consistency": 0.10,
            "RF": 0.10,
        }
    logger.info(
        "Weights (OD, BLV, GD, SO, GD_consistency, RF) = %.2f, %.2f, %.2f, %.2f, %.2f, %.2f",
        stage1_weights.get("deep", 0.0),
        stage1_weights.get("var", 0.0),
        stage1_weights.get("GD", 0.0),
        stage1_weights.get("species_overlap", 0.0),
        stage1_weights.get("gd_consistency", 0.0),
        stage1_weights.get("RF", 0.0),
    )


    try:
        for tree_id, tree_path in tree_dict.items():
            pbar.set_description(f"Processing {tree_id}")

            # 1. Read and initial processing
            Phylo_t0 = read_phylo_tree(tree_path)
            Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)
            Phylo_t2 = rename_input_tre(Phylo_t1, gene_to_new_name)

            # 2. Handle simple cases
            if len(get_species_set(Phylo_t2)) == 1 or len(get_species_list(Phylo_t2)) <= 3:
                rename_output_tre(
                    Phylo_t2,
                    new_name_to_gene,
                    tree_id,
                    dir_path,
                    source_tree_path=tree_path,
                )
                pbar.update(1)
                continue

            # 3. Depth annotation & Basal filtering
            Phylo_t2 = annotate_gene_tree(Phylo_t2, renamed_species_tree)

            current_gene_species = set(get_species_list(Phylo_t2))
            dynamic_basal_set = get_dynamic_basal_set(current_gene_species, renamed_species_tree)
            outgroup_root_list = get_all_rerooted_trees_filtered(
                Phylo_t2, dynamic_basal_set, new_name_to_gene
            )

            gd_nodes = find_dup_node(Phylo_t2, renamed_species_tree)
            gd_root_list = _reroot_by_node_set(Phylo_t2, gd_nodes)
            mad_root_list, minvar_root_list = _select_statistical_root_candidates(
                Phylo_t2, max_candidates=3
            )

            root_list = merge_unique_root_candidates(
                outgroup_root_list,
                gd_root_list,
                mad_root_list,
                minvar_root_list,
            )
            logger.debug(
                (
                    "%s: candidate roots from outgroup=%d, from GD_nodes=%d, "
                    "from MAD=%d, from MinVar=%d, merged=%d"
                ),
                tree_id,
                len(outgroup_root_list),
                len(gd_root_list),
                len(mad_root_list),
                len(minvar_root_list),
                len(root_list),
            )

            if not root_list:
                logger.warning("No valid root candidates after basal filter.")
                rename_output_tre(
                    Phylo_t2,
                    new_name_to_gene,
                    tree_id,
                    dir_path,
                    source_tree_path=tree_path,
                )
                pbar.update(1)
                continue

            # ==========================================
            # Single-stage scoring: compute all metrics including RF
            # ==========================================
            tree_objects = {}
            temp_stats = []

            for n, tree in enumerate(root_list):
                # tree is already in renamed Taxa ID format
                tree.resolve_polytomy(recursive=True)
                tree_key = f"{tree_id}_{n+1}"
                tree_objects[tree_key] = tree
                deep, var, GD, species_overlap, gd_consistency, RF = calculate_tree_statistics(
                    tree=tree,
                    species_tree=renamed_species_tree,
                )
                temp_stats.append({
                    "Tree": tree_key,
                    "deep": deep, "var": var, "GD": GD,
                    "species_overlap": species_overlap,
                    "gd_consistency": gd_consistency,
                    "RF": RF,
                    "tree_obj_ref": tree
                })
                # print(rename_input_tre(tree, new_name_to_gene))
                # print(deep, var, GD, species_overlap, gd_consistency, RF)
            current_df = pd.DataFrame(temp_stats)
            current_df["score"] = normalize_and_score(current_df, stage1_weights, include_rf=True)
            best_row = current_df.nlargest(1, "score").iloc[0]

            # ==========================================
            # Output
            # ==========================================
            best_tree_key = best_row["Tree"]
            best_tree = tree_objects[best_tree_key]

            rename_output_tre(
                best_tree,
                new_name_to_gene,
                tree_id,
                dir_path,
                source_tree_path=tree_path,
            )

            record = best_row.to_dict()
            del record["tree_obj_ref"]
            stat_matrix.append(record)

            pbar.update(1)

        # Save report
        stat_df = pd.DataFrame(stat_matrix)
        if not stat_df.empty:
            stat_df["tree_id"] = stat_df["Tree"].apply(lambda x: "_".join(x.split("_")[:-1]))
            cols = ["Tree", "score", "deep", "var", "GD", "species_overlap", "gd_consistency", "RF"]
            stat_df = stat_df.sort_values(by=["tree_id", "score"], ascending=[True, False])[cols]
            stat_df.to_csv(os.path.join(dir_path, "stat_matrix.csv"), index=False)  # save to dir_path not CWD

    finally:
        pbar.close()
