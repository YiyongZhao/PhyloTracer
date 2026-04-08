"""
MulRF distance calculation module for PhyloTracer.

This module computes species-level MulRF-like distances between multi-copy
gene trees and a reference species tree.
"""

import csv
import logging
from typing import Any, Dict, List, Optional, Set

from ete3 import Tree

from phylotracer import read_and_return_dict

logger = logging.getLogger(__name__)


def load_gene2sp_map(filepath: str) -> Dict[str, str]:
    """Load gene-to-species mapping from a two-column text file."""
    mapping = read_and_return_dict(filepath)
    return mapping


def infer_species_from_name(gene_name: str, separator: str = "_", position: str = "last") -> str:
    """Infer species from gene name using the selected separator rule."""
    if position == "last":
        parts = gene_name.rsplit(separator, 1)
        return parts[0] if len(parts) > 1 else gene_name
    parts = gene_name.split(separator, 1)
    return parts[0] if len(parts) > 1 else gene_name  # position=first should return prefix, not suffix


def build_species_map(
    tree: Tree,
    gene2sp_map: Optional[Dict[str, str]] = None,
    separator: str = "_",
    position: str = "last",
) -> Dict[str, str]:
    """Build {gene_leaf: species} mapping for a gene tree."""
    mapping = {}
    for leaf in tree.get_leaves():
        name = leaf.name
        if gene2sp_map and name in gene2sp_map:
            mapping[name] = gene2sp_map[name]
        else:
            mapping[name] = infer_species_from_name(name, separator, position)
    return mapping


def _precompute_gene_tree_features(
    gene_tree: Tree,
    gene2sp_map: Optional[Dict[str, str]] = None,
    separator: str = "_",
    position: str = "last",
) -> Dict[str, Any]:
    """Precompute reusable species-level structures for one gene tree.

    Returns a lightweight dict suitable for large-scale pairwise reuse.
    """
    sp_map = build_species_map(gene_tree, gene2sp_map, separator, position)
    species_set: Set[str] = set(sp_map.values())
    n_species = len(species_set)

    species_clades: Set[frozenset] = set()
    # Postorder dynamic aggregation avoids repeated node.get_leaves() calls.
    for node in gene_tree.traverse("postorder"):
        if node.is_leaf():
            node.add_feature("_sp_set_cache", {sp_map.get(node.name)})
            continue

        sp_union: Set[str] = set()
        for child in node.children:
            sp_union.update(getattr(child, "_sp_set_cache", set()))
        node.add_feature("_sp_set_cache", sp_union)

        if not node.is_root():
            clade = frozenset(sp for sp in sp_union if sp is not None)
            if 0 < len(clade) < n_species:
                species_clades.add(clade)

    return {
        "leaf_count": len(gene_tree.get_leaves()),
        "species_set": species_set,
        "species_count": n_species,
        "species_clades": species_clades,
    }


def _restrict_bipartitions_to_shared(
    species_clades: Set[frozenset],
    shared_species: Set[str],
) -> Set[frozenset]:
    """Restrict precomputed species clades to the current shared-species universe."""
    restricted: Set[frozenset] = set()
    n_shared = len(shared_species)
    for clade in species_clades:
        sub = clade & shared_species
        if 0 < len(sub) < n_shared:
            restricted.add(frozenset(sub))
    return restricted


def _empty_mode1_result(
    gt1_leaf_count: Optional[int] = None,
    gt2_leaf_count: Optional[int] = None,
    gt1_species_count: Optional[int] = None,
    gt2_species_count: Optional[int] = None,
    shared_species_count: Optional[int] = None,
) -> Dict[str, object]:
    """Build a mode=1-compatible empty result row."""
    return {
        "gene_tree_1_leaf_count": gt1_leaf_count,
        "gene_tree_2_leaf_count": gt2_leaf_count,
        "gene_tree_1_species_count": gt1_species_count,
        "gene_tree_2_species_count": gt2_species_count,
        "shared_species_count": shared_species_count,
        "mulrf_distance": None,
        "maximum_possible_mulrf_distance": None,
        "normalized_mulrf_distance": None,
        "shared_species_bipartition_count": None,
        "gene_tree_1_only_bipartition_count": None,
        "gene_tree_2_only_bipartition_count": None,
    }


def get_bipartitions_species_level(gene_tree: Tree, sp_map: Dict[str, str], shared_species: set) -> set:
    """Extract species-level bipartitions from a multi-copy gene tree."""
    bipartitions = set()
    all_sp = frozenset(shared_species)
    for node in gene_tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        clade_species = frozenset(
            sp_map[l.name]
            for l in node.get_leaves()
            if sp_map.get(l.name) in shared_species
        )
        if 0 < len(clade_species) < len(all_sp):
            bipartitions.add(clade_species)
    return bipartitions


def get_bipartitions_species_tree(sp_tree: Tree, shared_species: set) -> set:
    """Extract bipartitions from a species tree."""
    bipartitions = set()
    all_sp = frozenset(shared_species)
    for node in sp_tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        clade_species = frozenset(
            l.name for l in node.get_leaves() if l.name in shared_species
        )
        if 0 < len(clade_species) < len(all_sp):
            bipartitions.add(clade_species)
    return bipartitions


def compute_mulrf(
    gene_tree: Tree,
    sp_tree: Tree,
    gene2sp_map: Optional[Dict[str, str]] = None,
    separator: str = "_",
    position: str = "last",
) -> Dict[str, object]:
    """Compute MulRF distance for one gene tree."""
    sp_map = build_species_map(gene_tree, gene2sp_map, separator, position)

    gene_species = set(sp_map.values())
    sp_tree_species = set(sp_tree.get_leaf_names())
    shared_species = gene_species & sp_tree_species

    n_leaves = len(gene_tree.get_leaves())
    n_species = len(gene_species)
    n_shared = len(shared_species)

    if n_shared < 2:
        return {
            "n_leaves": n_leaves,
            "n_species": n_species,
            "n_shared_species": n_shared,
            "mulrf": None,
            "max_rf": None,
            "normalized_mulrf": None,
            "shared_bipartitions": None,
            "only_in_gene": None,
            "only_in_sp": None,
            "error": f"Only {n_shared} shared species (need >=2)",
        }

    pruned_sp = sp_tree.copy("newick")
    keep = sorted(shared_species & set(pruned_sp.get_leaf_names()))
    if len(keep) < 2:
        return {
            "n_leaves": n_leaves,
            "n_species": n_species,
            "n_shared_species": n_shared,
            "mulrf": None,
            "max_rf": None,
            "normalized_mulrf": None,
            "shared_bipartitions": None,
            "only_in_gene": None,
            "only_in_sp": None,
            "error": "Too few shared species after pruning",
        }
    if len(keep) != len(pruned_sp.get_leaf_names()):
        pruned_sp.prune(keep, preserve_branch_length=True)

    gene_bip = get_bipartitions_species_level(gene_tree, sp_map, shared_species)
    sp_bip = get_bipartitions_species_tree(pruned_sp, shared_species)

    shared = gene_bip & sp_bip
    only_gene = gene_bip - sp_bip
    only_sp = sp_bip - gene_bip

    mulrf = len(only_gene) + len(only_sp)
    max_rf = len(gene_bip) + len(sp_bip)
    normalized = mulrf / max_rf if max_rf > 0 else 0.0

    return {
        "n_leaves": n_leaves,
        "n_species": n_species,
        "n_shared_species": n_shared,
        "mulrf": mulrf,
        "max_rf": max_rf,
        "normalized_mulrf": round(normalized, 6),
        "shared_bipartitions": len(shared),
        "only_in_gene": len(only_gene),
        "only_in_sp": len(only_sp),
        "error": "",
    }


def compute_mulrf_gene_vs_species(
    gene_tree: Tree,
    sp_tree: Tree,
    gene2sp_map: Optional[Dict[str, str]] = None,
    separator: str = "_",
    position: str = "last",
) -> Dict[str, object]:
    """Compute MulRF distance between a gene tree and a species tree.

    This is a convenience wrapper around :func:`compute_mulrf` that returns
    results using keys consistent with :func:`compute_mulrf_between_gene_trees`.
    """
    res = compute_mulrf(gene_tree, sp_tree, gene2sp_map, separator, position)
    mulrf_dist = res.get("mulrf")
    return {
        "gene_tree_leaf_count": res.get("n_leaves"),
        "gene_tree_species_count": res.get("n_species"),
        "shared_species_count": res.get("n_shared_species"),
        "mulrf_distance": mulrf_dist,
        "maximum_possible_mulrf_distance": res.get("max_rf"),
        "normalized_mulrf_distance": res.get("normalized_mulrf"),
        "shared_bipartitions": res.get("shared_bipartitions"),
        "only_in_gene": res.get("only_in_gene"),
        "only_in_sp": res.get("only_in_sp"),
        "error": res.get("error", ""),
    }


def compute_mulrf_between_gene_trees(
    gene_tree_1: Optional[Tree],
    gene_tree_2: Optional[Tree],
    gene2sp_map_1: Optional[Dict[str, str]] = None,
    gene2sp_map_2: Optional[Dict[str, str]] = None,
    separator: str = "_",
    position: str = "last",
    precomputed_1: Optional[Dict[str, Any]] = None,
    precomputed_2: Optional[Dict[str, Any]] = None,
) -> Dict[str, object]:
    """Compute species-level MulRF distance between two gene trees.

    If precomputed_* dicts are provided, no tree parsing/traversal is repeated.
    """
    if precomputed_1 is None:
        if gene_tree_1 is None:
            raise ValueError("gene_tree_1 is required when precomputed_1 is not provided")
        precomputed_1 = _precompute_gene_tree_features(gene_tree_1, gene2sp_map_1, separator, position)

    if precomputed_2 is None:
        if gene_tree_2 is None:
            raise ValueError("gene_tree_2 is required when precomputed_2 is not provided")
        precomputed_2 = _precompute_gene_tree_features(gene_tree_2, gene2sp_map_2, separator, position)

    species_1 = precomputed_1["species_set"]
    species_2 = precomputed_2["species_set"]
    shared_species = species_1 & species_2

    if len(shared_species) < 2:
        return _empty_mode1_result(
            gt1_leaf_count=precomputed_1["leaf_count"],
            gt2_leaf_count=precomputed_2["leaf_count"],
            gt1_species_count=precomputed_1["species_count"],
            gt2_species_count=precomputed_2["species_count"],
            shared_species_count=len(shared_species),
        )

    bip_1 = _restrict_bipartitions_to_shared(precomputed_1["species_clades"], shared_species)
    bip_2 = _restrict_bipartitions_to_shared(precomputed_2["species_clades"], shared_species)

    shared = bip_1 & bip_2
    only_1 = bip_1 - bip_2
    only_2 = bip_2 - bip_1
    mulrf = len(only_1) + len(only_2)
    max_rf = len(bip_1) + len(bip_2)
    normalized = mulrf / max_rf if max_rf > 0 else 0.0

    return {
        "gene_tree_1_leaf_count": precomputed_1["leaf_count"],
        "gene_tree_2_leaf_count": precomputed_2["leaf_count"],
        "gene_tree_1_species_count": precomputed_1["species_count"],
        "gene_tree_2_species_count": precomputed_2["species_count"],
        "shared_species_count": len(shared_species),
        "mulrf_distance": mulrf,
        "maximum_possible_mulrf_distance": max_rf,
        "normalized_mulrf_distance": round(normalized, 6),
        "shared_species_bipartition_count": len(shared),
        "gene_tree_1_only_bipartition_count": len(only_1),
        "gene_tree_2_only_bipartition_count": len(only_2),
    }


def mulrf_main(
    tre_dic: Dict[str, str],
    species_tree_path: Optional[str] = None,
    gene2sp_map_path: Optional[str] = None,
    mode: str = "1",
    output_file: str = "mulrf_distance.tsv",
    separator: str = "_",
    position: str = "last",
) -> None:
    """
    Batch MulRF calculation.

    mode = "1": Gene Tree vs Gene Tree (pairwise within one GF list, excluding self-comparisons)
    mode = "2": Gene Tree vs Species Tree
    """
    gene2sp_map_1 = load_gene2sp_map(gene2sp_map_path) if gene2sp_map_path else None

    if mode == "1":
        columns = [
            "tre_id_1",
            "tre_id_2",
            "gene_tree_1_leaf_count",
            "gene_tree_2_leaf_count",
            "gene_tree_1_species_count",
            "gene_tree_2_species_count",
            "shared_species_count",
            "mulrf_distance",
            "maximum_possible_mulrf_distance",
            "normalized_mulrf_distance",
            "shared_species_bipartition_count",
            "gene_tree_1_only_bipartition_count",
            "gene_tree_2_only_bipartition_count",
        ]

        items = list(tre_dic.items())
        n_trees = len(items)
        total_pairs = n_trees * (n_trees - 1) // 2
        logger.info("MulRF mode=1: total_trees=%d", n_trees)
        logger.info("MulRF mode=1: total_expected_pairs=%d", total_pairs)

        precomputed: Dict[str, Optional[Dict[str, Any]]] = {}
        for idx, (tre_id, tree_path) in enumerate(items, start=1):
            try:
                gt = Tree(tree_path)
                precomputed[tre_id] = _precompute_gene_tree_features(
                    gt,
                    gene2sp_map_1,
                    separator,
                    position,
                )
            except Exception as e:
                logger.warning("Failed to parse tree %s (%s): %s", tre_id, tree_path, e)
                precomputed[tre_id] = None

            if idx % 1000 == 0:
                logger.info("Precompute progress: %d/%d trees", idx, n_trees)

        progress_interval = 100000
        processed_pairs = 0
        valid_pairs = 0
        mulrf_sum = 0.0

        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
            writer.writeheader()

            for i in range(n_trees):
                tre_id_1, _ = items[i]
                pre1 = precomputed.get(tre_id_1)

                for j in range(i + 1, n_trees):
                    tre_id_2, _ = items[j]
                    pre2 = precomputed.get(tre_id_2)

                    if pre1 is None or pre2 is None:
                        res = _empty_mode1_result(
                            gt1_leaf_count=(pre1 or {}).get("leaf_count"),
                            gt2_leaf_count=(pre2 or {}).get("leaf_count"),
                            gt1_species_count=(pre1 or {}).get("species_count"),
                            gt2_species_count=(pre2 or {}).get("species_count"),
                            shared_species_count=None,
                        )
                    else:
                        try:
                            res = compute_mulrf_between_gene_trees(
                                None,
                                None,
                                gene2sp_map_1,
                                gene2sp_map_1,
                                separator,
                                position,
                                precomputed_1=pre1,
                                precomputed_2=pre2,
                            )
                        except Exception as e:
                            logger.warning(
                                "MulRF computation failed for pair (%s, %s): %s",
                                tre_id_1,
                                tre_id_2,
                                e,
                            )
                            res = _empty_mode1_result(
                                gt1_leaf_count=pre1.get("leaf_count"),
                                gt2_leaf_count=pre2.get("leaf_count"),
                                gt1_species_count=pre1.get("species_count"),
                                gt2_species_count=pre2.get("species_count"),
                                shared_species_count=None,
                            )

                    row = {"tre_id_1": tre_id_1, "tre_id_2": tre_id_2}
                    row.update(res)
                    writer.writerow(row)

                    processed_pairs += 1
                    if row.get("mulrf_distance") is not None:
                        valid_pairs += 1
                        mulrf_sum += row["mulrf_distance"]

                    if processed_pairs % progress_interval == 0:
                        logger.info(
                            "MulRF mode=1 progress: %d/%d pairs (%.2f%%)",
                            processed_pairs,
                            total_pairs,
                            (processed_pairs / total_pairs * 100.0) if total_pairs > 0 else 100.0,
                        )

        mean_mulrf = mulrf_sum / valid_pairs if valid_pairs else 0.0
        logger.info("MulRF pairwise summary: valid_pairs=%d, mean_mulrf=%.4f", valid_pairs, mean_mulrf)
        logger.info("MulRF results written to %s", output_file)
        return

    rows: List[Dict[str, object]] = []

    if mode == "2":
        if not species_tree_path:
            raise ValueError("species_tree_path is required when mode=2.")
        sp_tree = Tree(species_tree_path)
        for tre_id, tree_path in tre_dic.items():
            try:
                gt = Tree(tree_path)
                res = compute_mulrf(gt, sp_tree, gene2sp_map_1, separator, position)
                row = {
                    "tre_id": tre_id,
                    "gene_tree_leaf_count": res.get("n_leaves"),
                    "gene_tree_species_count": res.get("n_species"),
                    "mulrf_distance": res.get("mulrf"),
                    "maximum_possible_mulrf_distance": res.get("max_rf"),
                    "normalized_mulrf_distance": res.get("normalized_mulrf"),
                    "shared_species_bipartition_count": res.get("shared_bipartitions"),
                    "gene_tree_only_bipartition_count": res.get("only_in_gene"),
                }
            except Exception as e:
                logger.warning("MulRF computation failed for tree pair: %s", e)
                row = {
                    "tre_id": tre_id,
                    "gene_tree_leaf_count": None,
                    "gene_tree_species_count": None,
                    "mulrf_distance": None,
                    "maximum_possible_mulrf_distance": None,
                    "normalized_mulrf_distance": None,
                    "shared_species_bipartition_count": None,
                    "gene_tree_only_bipartition_count": None,
                }
            rows.append(row)
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    rows.sort(key=lambda x: (x["mulrf_distance"] if x["mulrf_distance"] is not None else float("inf")))
    columns = [
        "tre_id",
        "gene_tree_leaf_count",
        "gene_tree_species_count",
        "mulrf_distance",
        "maximum_possible_mulrf_distance",
        "normalized_mulrf_distance",
        "shared_species_bipartition_count",
        "gene_tree_only_bipartition_count",
    ]
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    valid = [r for r in rows if r.get("mulrf_distance") is not None]
    mean_mulrf = sum(r["mulrf_distance"] for r in valid) / len(valid) if valid else 0.0
    logger.info("MulRF pairwise summary: valid_pairs=%d, mean_mulrf=%.4f", len(valid), mean_mulrf)
    logger.info("MulRF results written to %s", output_file)
