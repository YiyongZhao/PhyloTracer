#!/usr/bin/env python3
"""
Benchmark script for PhyloTracer bioinformatics toolkit.

Measures wall-clock time and peak memory usage for key PhyloTracer operations
using synthetic test data (random trees generated via ete3). No external data
files are required.

Usage:
    python benchmark.py
    python benchmark.py --sizes 50 100 500
    python benchmark.py --sizes 50 100 500 1000 --repeat 3
"""

import argparse
import gc
import random
import sys
import time
import tracemalloc
from contextlib import contextmanager

# ---------------------------------------------------------------------------
# Graceful import checking
# ---------------------------------------------------------------------------

_MISSING_DEPS = []

try:
    from ete3 import Tree, PhyloTree
except ImportError:
    _MISSING_DEPS.append("ete3")

try:
    import numpy as np
except ImportError:
    _MISSING_DEPS.append("numpy")

if _MISSING_DEPS:
    print(
        "ERROR: The following required packages are not installed: "
        + ", ".join(_MISSING_DEPS),
        file=sys.stderr,
    )
    print(
        "Install them with:  pip install " + " ".join(_MISSING_DEPS),
        file=sys.stderr,
    )
    sys.exit(1)

# Try importing PhyloTracer helpers; fall back to local implementations if the
# package is not installed.
try:
    from phylotracer import (
        get_species_set,
        get_species_list,
        annotate_gene_tree,
        find_dup_node,
        root_tre_with_midpoint_outgroup,
        is_rooted,
    )
    _PHYLOTRACER_AVAILABLE = True
except ImportError:
    _PHYLOTRACER_AVAILABLE = False


# ===================================================================
# Fallback implementations (used only when PhyloTracer is not installed)
# ===================================================================

if not _PHYLOTRACER_AVAILABLE:
    def get_species_list(node):
        if node is None:
            return []
        return [leaf.name.split("_")[0] for leaf in node.iter_leaves()]

    def get_species_set(node):
        return set(get_species_list(node))

    def is_rooted(tree):
        return len(tree.get_children()) == 2

    def root_tre_with_midpoint_outgroup(tree):
        tree_copy = tree.copy("newick")
        if is_rooted(tree_copy):
            return tree_copy
        leaves = tree_copy.get_leaves()
        if len(leaves) < 3:
            return tree_copy
        try:
            mid = tree_copy.get_midpoint_outgroup()
            tree_copy.set_outgroup(mid if not mid.is_root() else leaves[0])
        except Exception:
            tree_copy.set_outgroup(leaves[0])
        return tree_copy

    def annotate_gene_tree(gene_tree, species_tree):
        for s_node in species_tree.traverse("preorder"):
            if s_node.up is None:
                s_node.add_feature("depth", 0)
            else:
                s_node.add_feature("depth", s_node.up.depth + 1)
        for node in gene_tree.traverse("postorder"):
            curr_sps = get_species_set(node)
            if not curr_sps:
                continue
            try:
                if len(curr_sps) == 1:
                    mapped_node = species_tree & list(curr_sps)[0]
                else:
                    mapped_node = species_tree.get_common_ancestor(list(curr_sps))
                node.add_feature("map", mapped_node.name)
                node.add_feature("depth", mapped_node.depth)
            except Exception:
                node.add_feature("map", "unknown")
                node.add_feature("depth", None)
            if not node.is_leaf():
                children = node.get_children()
                if len(children) == 2:
                    sps_a = get_species_set(children[0])
                    sps_b = get_species_set(children[1])
                    overlap = len(sps_a & sps_b)
                    node.add_feature("overlap", overlap)
                    node.add_feature("is_gd", overlap >= 1)
                else:
                    node.add_feature("overlap", 0)
                    node.add_feature("is_gd", False)
            else:
                node.add_feature("overlap", 0)
                node.add_feature("is_gd", False)
        return gene_tree

    def find_dup_node(gene_tree, species_tree, gd_support=50, clade_support=50,
                      dup_species_num=1, dup_species_percent=0,
                      max_topology_distance=1):
        dup_nodes = []
        for node in gene_tree.traverse("postorder"):
            if node.is_leaf() or not getattr(node, "is_gd", False):
                continue
            children = node.get_children()
            if len(children) != 2:
                continue
            sps_a = get_species_set(children[0])
            sps_b = get_species_set(children[1])
            overlap_sps = sps_a & sps_b
            if len(overlap_sps) >= dup_species_num:
                dup_nodes.append(node)
        return dup_nodes


# ===================================================================
# Synthetic data generation
# ===================================================================

SPECIES_NAMES = [
    "Arabidopsis", "Oryza", "Zea", "Solanum", "Vitis",
    "Populus", "Glycine", "Medicago", "Brassica", "Gossypium",
    "Malus", "Prunus", "Citrus", "Eucalyptus", "Theobroma",
    "Coffea", "Helianthus", "Lactuca", "Daucus", "Beta",
    "Spinacia", "Cucumis", "Citrullus", "Capsicum", "Nicotiana",
    "Petunia", "Ipomoea", "Manihot", "Ricinus", "Jatropha",
    "Hevea", "Linum", "Cannabis", "Humulus", "Morus",
    "Ficus", "Juglans", "Quercus", "Fagus", "Betula",
]


def _generate_species_tree(n_species):
    """Create a random species tree with the given number of taxa."""
    sptree = Tree()
    sptree.populate(n_species)
    species = SPECIES_NAMES[:n_species] if n_species <= len(SPECIES_NAMES) else [
        f"Sp{i}" for i in range(n_species)
    ]
    for i, leaf in enumerate(sptree.get_leaves()):
        leaf.name = species[i]
    # Label internal nodes
    idx = 0
    for node in sptree.traverse("postorder"):
        if not node.is_leaf():
            node.name = f"S{idx}"
            idx += 1
    return sptree, species


def _generate_gene_tree(n_taxa, species_list):
    """Create a random gene tree whose leaves are named Species_GeneN."""
    tree = Tree()
    tree.populate(n_taxa, random_branches=True)
    n_species = len(species_list)
    for i, leaf in enumerate(tree.get_leaves()):
        sp = species_list[i % n_species]
        gene_num = (i // n_species) + 1
        leaf.name = f"{sp}_{gene_num}"
        # Assign a random support value to the parent if it exists
        if leaf.up is not None and leaf.up.support == 1.0:
            leaf.up.support = random.uniform(30, 100)
    # Set support values on internal nodes
    for node in tree.traverse():
        if not node.is_leaf() and not node.is_root():
            node.support = random.uniform(30, 100)
    return tree


def _generate_newick_string(n_taxa):
    """Return a Newick string for a random tree with n_taxa leaves."""
    tree = Tree()
    tree.populate(n_taxa, random_branches=True)
    return tree.write(format=5)


# ===================================================================
# Measurement utilities
# ===================================================================

@contextmanager
def measure():
    """Context manager that yields a dict to be filled with time_s and peak_mb."""
    gc.collect()
    tracemalloc.start()
    start_time = time.perf_counter()
    result = {}
    try:
        yield result
    finally:
        elapsed = time.perf_counter() - start_time
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        result["time_s"] = elapsed
        result["peak_mb"] = peak / (1024 * 1024)


# ===================================================================
# Benchmark functions
# ===================================================================

def bench_tree_read_parse(n_taxa, repeat):
    """Benchmark reading and parsing a Newick tree string."""
    newick = _generate_newick_string(n_taxa)
    times = []
    peaks = []
    for _ in range(repeat):
        with measure() as m:
            _ = Tree(newick, format=5)
        times.append(m["time_s"])
        peaks.append(m["peak_mb"])
    return {
        "operation": "Tree read/parse",
        "input_size": n_taxa,
        "time_s": sum(times) / len(times),
        "peak_mb": max(peaks),
    }


def bench_tree_rooting(n_taxa, repeat):
    """Benchmark midpoint rooting of a tree."""
    species_list = SPECIES_NAMES[:min(n_taxa, len(SPECIES_NAMES))]
    if n_taxa > len(SPECIES_NAMES):
        species_list = [f"Sp{i}" for i in range(n_taxa)]
    tree = _generate_gene_tree(n_taxa, species_list)
    times = []
    peaks = []
    for _ in range(repeat):
        with measure() as m:
            _ = root_tre_with_midpoint_outgroup(tree)
        times.append(m["time_s"])
        peaks.append(m["peak_mb"])
    return {
        "operation": "Tree rooting (midpoint)",
        "input_size": n_taxa,
        "time_s": sum(times) / len(times),
        "peak_mb": max(peaks),
    }


def bench_species_set_extraction(n_taxa, repeat):
    """Benchmark extracting the species set from a gene tree."""
    n_species = min(n_taxa // 2, len(SPECIES_NAMES)) or 5
    species_list = SPECIES_NAMES[:n_species]
    tree = _generate_gene_tree(n_taxa, species_list)
    times = []
    peaks = []
    for _ in range(repeat):
        with measure() as m:
            for node in tree.traverse():
                _ = get_species_set(node)
        times.append(m["time_s"])
        peaks.append(m["peak_mb"])
    return {
        "operation": "Species set extraction",
        "input_size": n_taxa,
        "time_s": sum(times) / len(times),
        "peak_mb": max(peaks),
    }


def bench_gene_duplication_detection(n_taxa, repeat):
    """Benchmark gene duplication detection via annotation and find_dup_node."""
    n_species = min(n_taxa // 2, len(SPECIES_NAMES)) or 5
    species_list = SPECIES_NAMES[:n_species]
    sptree, _ = _generate_species_tree(n_species)
    gene_tree = _generate_gene_tree(n_taxa, species_list)
    # Ensure the gene tree is rooted for annotation
    gene_tree = root_tre_with_midpoint_outgroup(gene_tree)
    times = []
    peaks = []
    for _ in range(repeat):
        gt_copy = gene_tree.copy("newick")
        # Re-parse as PhyloTree for compatibility
        gt_copy = PhyloTree(gt_copy.write(format=5), format=5)
        sp_copy = sptree.copy("newick")
        with measure() as m:
            annotate_gene_tree(gt_copy, sp_copy)
            _ = find_dup_node(gt_copy, sp_copy)
        times.append(m["time_s"])
        peaks.append(m["peak_mb"])
    return {
        "operation": "Gene duplication detection",
        "input_size": n_taxa,
        "time_s": sum(times) / len(times),
        "peak_mb": max(peaks),
    }


def bench_tree_topology_comparison(n_taxa, repeat):
    """Benchmark Robinson-Foulds topology comparison between two random trees."""
    n_species = min(n_taxa, len(SPECIES_NAMES)) or 5
    species_list = SPECIES_NAMES[:n_species]
    if n_taxa > len(SPECIES_NAMES):
        species_list = [f"Sp{i}" for i in range(n_taxa)]

    tree_a = Tree()
    tree_a.populate(n_taxa)
    tree_b = Tree()
    tree_b.populate(n_taxa)
    # Ensure both trees share the same leaf names
    leaf_names = [f"{species_list[i % len(species_list)]}_{i}" for i in range(n_taxa)]
    for i, leaf in enumerate(tree_a.get_leaves()):
        leaf.name = leaf_names[i]
    for i, leaf in enumerate(tree_b.get_leaves()):
        leaf.name = leaf_names[i]

    times = []
    peaks = []
    for _ in range(repeat):
        with measure() as m:
            _ = tree_a.robinson_foulds(tree_b, unrooted_trees=True)
        times.append(m["time_s"])
        peaks.append(m["peak_mb"])
    return {
        "operation": "Tree topology comparison (RF)",
        "input_size": n_taxa,
        "time_s": sum(times) / len(times),
        "peak_mb": max(peaks),
    }


# ===================================================================
# Output formatting
# ===================================================================

def print_results_table(results):
    """Print a formatted table of benchmark results."""
    header = f"{'Operation':<35} {'Input Size':>10} {'Time (s)':>12} {'Peak Mem (MB)':>14}"
    separator = "-" * len(header)

    print()
    print("=" * len(header))
    print("PhyloTracer Benchmark Results")
    print("=" * len(header))
    if _PHYLOTRACER_AVAILABLE:
        print("(using installed PhyloTracer package)")
    else:
        print("(using built-in fallback implementations)")
    print()
    print(header)
    print(separator)

    current_op = None
    for r in results:
        if r["operation"] != current_op:
            if current_op is not None:
                print(separator)
            current_op = r["operation"]
        print(
            f"{r['operation']:<35} {r['input_size']:>10} "
            f"{r['time_s']:>12.6f} {r['peak_mb']:>14.3f}"
        )

    print(separator)
    print()


# ===================================================================
# Main entry point
# ===================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark key PhyloTracer operations using synthetic tree data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python benchmark.py\n"
            "  python benchmark.py --sizes 50 100 500\n"
            "  python benchmark.py --sizes 50 100 500 1000 --repeat 5\n"
        ),
    )
    parser.add_argument(
        "--sizes",
        type=int,
        nargs="+",
        default=[50, 100, 500, 1000],
        help="Tree sizes (number of taxa) to benchmark (default: 50 100 500 1000)",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=3,
        help="Number of repetitions per measurement for averaging (default: 3)",
    )
    args = parser.parse_args()

    sizes = sorted(args.sizes)
    repeat = args.repeat

    print(f"Running benchmarks with sizes={sizes}, repeat={repeat}")
    print(f"Python {sys.version}")

    benchmarks = [
        ("Tree read/parse", bench_tree_read_parse),
        ("Tree rooting (midpoint)", bench_tree_rooting),
        ("Species set extraction", bench_species_set_extraction),
        ("Gene duplication detection", bench_gene_duplication_detection),
        ("Tree topology comparison (RF)", bench_tree_topology_comparison),
    ]

    all_results = []

    for bname, bfunc in benchmarks:
        for size in sizes:
            sys.stdout.write(f"  {bname} (n={size}) ... ")
            sys.stdout.flush()
            try:
                result = bfunc(size, repeat)
                all_results.append(result)
                sys.stdout.write(f"done ({result['time_s']:.4f}s)\n")
            except Exception as exc:
                sys.stdout.write(f"FAILED: {exc}\n")
                all_results.append({
                    "operation": bname,
                    "input_size": size,
                    "time_s": float("nan"),
                    "peak_mb": float("nan"),
                })

    print_results_table(all_results)


if __name__ == "__main__":
    main()
