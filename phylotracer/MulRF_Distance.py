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


def mulrf_main(
    tre_dic: Dict[str, str],
    species_tree_path: str,
    gene2sp_map_path: Optional[str] = None,
    separator: str = "_",
    position: str = "last",
    quiet: bool = False,
) -> None:
    """Batch MulRF calculation for all gene trees in ``tre_dic``."""
    gene2sp_map = load_gene2sp_map(gene2sp_map_path) if gene2sp_map_path else None
    sp_tree = Tree(species_tree_path)

    rows = []
    for idx, (family, tree_path) in enumerate(tre_dic.items()):
        try:
            gt = Tree(tree_path)
            res = compute_mulrf(gt, sp_tree, gene2sp_map, separator, position)
        except Exception as exc:
            res = {
                "n_leaves": None,
                "n_species": None,
                "n_shared_species": None,
                "mulrf": None,
                "max_rf": None,
                "normalized_mulrf": None,
                "shared_bipartitions": None,
                "only_in_gene": None,
                "only_in_sp": None,
                "error": f"Parse/compute error: {exc}",
            }
        res["tree_index"] = idx
        res["gene_family"] = family
        rows.append(res)

    rows.sort(key=lambda x: (x["mulrf"] if x["mulrf"] is not None else float("inf")))

    out_file = "mulrf_distance.tsv"
    columns = [
        "tree_index",
        "gene_family",
        "n_leaves",
        "n_species",
        "n_shared_species",
        "mulrf",
        "max_rf",
        "normalized_mulrf",
        "shared_bipartitions",
        "only_in_gene",
        "only_in_sp",
        "error",
    ]
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    valid = [r for r in rows if r.get("mulrf") is not None]
    if not quiet and valid:
        mean_mulrf = sum(r["mulrf"] for r in valid) / len(valid)
        logger.info("MulRF summary: valid_trees=%d, mean_mulrf=%.4f", len(valid), mean_mulrf)
    logger.info("MulRF results written to %s", out_file)
