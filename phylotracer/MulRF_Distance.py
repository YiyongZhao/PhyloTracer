"""
MulRF distance calculation module for PhyloTracer.

This module computes species-level MulRF-like distances between multi-copy
gene trees and a reference species tree.
"""

import csv
import logging
from typing import Dict, Optional

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
    return parts[1] if len(parts) > 1 else gene_name


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


def compute_mulrf_between_gene_trees(
    gene_tree_1: Tree,
    gene_tree_2: Tree,
    gene2sp_map_1: Optional[Dict[str, str]] = None,
    gene2sp_map_2: Optional[Dict[str, str]] = None,
    separator: str = "_",
    position: str = "last",
) -> Dict[str, object]:
    """Compute species-level MulRF distance between two gene trees."""
    sp_map_1 = build_species_map(gene_tree_1, gene2sp_map_1, separator, position)
    sp_map_2 = build_species_map(gene_tree_2, gene2sp_map_2, separator, position)

    species_1 = set(sp_map_1.values())
    species_2 = set(sp_map_2.values())
    shared_species = species_1 & species_2

    if len(shared_species) < 2:
        return {
            "gene_tree_1_leaf_count": len(gene_tree_1.get_leaves()),
            "gene_tree_2_leaf_count": len(gene_tree_2.get_leaves()),
            "gene_tree_1_species_count": len(species_1),
            "gene_tree_2_species_count": len(species_2),
            "shared_species_count": len(shared_species),
            "mulrf_distance": None,
            "maximum_possible_mulrf_distance": None,
            "normalized_mulrf_distance": None,
            "shared_species_bipartition_count": None,
            "gene_tree_1_only_bipartition_count": None,
            "gene_tree_2_only_bipartition_count": None,
        }

    bip_1 = get_bipartitions_species_level(gene_tree_1, sp_map_1, shared_species)
    bip_2 = get_bipartitions_species_level(gene_tree_2, sp_map_2, shared_species)

    shared = bip_1 & bip_2
    only_1 = bip_1 - bip_2
    only_2 = bip_2 - bip_1
    mulrf = len(only_1) + len(only_2)
    max_rf = len(bip_1) + len(bip_2)
    normalized = mulrf / max_rf if max_rf > 0 else 0.0

    return {
        "gene_tree_1_leaf_count": len(gene_tree_1.get_leaves()),
        "gene_tree_2_leaf_count": len(gene_tree_2.get_leaves()),
        "gene_tree_1_species_count": len(species_1),
        "gene_tree_2_species_count": len(species_2),
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

    mode = \"1\": Gene Tree vs Gene Tree (pairwise within one GF list, excluding self-comparisons)
    mode = \"2\": Gene Tree vs Species Tree
    """
    gene2sp_map_1 = load_gene2sp_map(gene2sp_map_path) if gene2sp_map_path else None
    rows = []

    if mode == "1":
        items = list(tre_dic.items())
        for i in range(len(items)):
            tre_id_1, tree_path_1 = items[i]
            for j in range(i + 1, len(items)):
                tre_id_2, tree_path_2 = items[j]
                try:
                    gt1 = Tree(tree_path_1)
                    gt2 = Tree(tree_path_2)
                    res = compute_mulrf_between_gene_trees(
                        gt1,
                        gt2,
                        gene2sp_map_1,
                        gene2sp_map_1,
                        separator,
                        position,
                    )
                except Exception as e:
                    logger.warning(f"MulRF computation failed for tree pair: {e}")
                    res = {
                        "gene_tree_1_leaf_count": None,
                        "gene_tree_2_leaf_count": None,
                        "gene_tree_1_species_count": None,
                        "gene_tree_2_species_count": None,
                        "shared_species_count": None,
                        "mulrf_distance": None,
                        "maximum_possible_mulrf_distance": None,
                        "normalized_mulrf_distance": None,
                        "shared_species_bipartition_count": None,
                        "gene_tree_1_only_bipartition_count": None,
                        "gene_tree_2_only_bipartition_count": None,
                    }
                row = {"tre_id_1": tre_id_1, "tre_id_2": tre_id_2}
                row.update(res)
                rows.append(row)
    elif mode == "2":
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
                logger.warning(f"MulRF computation failed for tree pair: {e}")
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
    if mode == "2":
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
    else:
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
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    valid = [r for r in rows if r.get("mulrf_distance") is not None]
    mean_mulrf = sum(r["mulrf_distance"] for r in valid) / len(valid) if valid else 0.0
    logger.info("MulRF pairwise summary: valid_pairs=%d, mean_mulrf=%.4f", len(valid), mean_mulrf)
    logger.info("MulRF results written to %s", output_file)
