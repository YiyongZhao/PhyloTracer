"""
Shared utilities and common imports for the PhyloTracer phylogenomics pipeline.

This module centralizes mapping, tree-processing, and duplication-detection
helpers used across the project to improve reproducibility and traceability.
"""

import logging
logger = logging.getLogger(__name__)
import os
import random
import shutil
from typing import Dict, List, Optional, Set, Tuple

import re
import numpy as np
import pandas as pd
from ete3 import PhyloTree, Tree
try:
    from ete3 import NodeStyle, TextFace, TreeStyle
except ImportError:
    NodeStyle = None
    TextFace = None
    TreeStyle = None
from tqdm import tqdm

__all__ = [
    "read_and_return_dict",
    "generate_sps_voucher",
    "gene_id_transfer",
    "rename_input_tre",
    "read_tree",
    "read_phylo_tree",
    "is_rooted",
    "root_tre_with_midpoint_outgroup",
    "num_sptree",
    "num_tre_node",
    "get_max_deepth",
    "calculate_depth",
    "calculate_deepvar",
    "get_species_list",
    "get_species_set",
    "calculate_species_num",
    "annotate_gene_tree",
    "compute_tip_to_root_branch_length_variance",
    "realign_branch_length",
    "rejust_root_dist",
    "judge_support",
    "sps_dup_num",
    "get_gene_pairs",
    "find_tre_dup",
    "map_species_set_to_node",
    "find_dup_node",
    "calculate_gd_num",
]

# =========================
# I/O & Mapping Utilities
# =========================


def read_and_return_dict(filename: str, separator: str = "\t") -> Dict[str, str]:
    """Parse a two-column mapping file into a dictionary.

    Args:
        filename (str): Path to a delimiter-separated file containing key/value
            pairs in the first two columns.
        separator (str): Column delimiter used in the file.

    Returns:
        Dict[str, str]: Mapping from column-1 keys to column-2 values.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is empty.
        Exception: If parsing fails for any other reason.

    Assumptions:
        The file has no header and at least two columns per row.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    try:
        df = pd.read_csv(filename, sep=separator, header=None)
        if df.empty:
            raise ValueError("Mapping file is empty")
        if 1 not in df.columns and separator == "\t":
            # Fallback for legacy mapping files using generic whitespace
            # delimiters instead of tabs.
            df = pd.read_csv(filename, sep=r"\s+", header=None, engine="python")
        if 1 not in df.columns:
            raise ValueError(
                "Mapping file must contain at least two columns "
                "(key and value)."
            )
        return df.set_index(0)[1].to_dict()
    except Exception as exc:
        logging.error(f"Failed to parse mapping file {filename}: {exc}")
        raise



def generate_sps_voucher(sps_num: int) -> List[str]:
    """Generate unique three-character vouchers for species labeling.

    Args:
        sps_num (int): Number of species vouchers required.

    Returns:
        List[str]: Sorted list of unique voucher strings.

    Assumptions:
        Vouchers are sampled without replacement from 62 alphanumeric
        characters; the maximum available unique vouchers is 62P3.
    """
    characters = (
        [chr(i) for i in range(65, 91)]
        + [chr(i) for i in range(97, 123)]
        + [str(i) for i in range(10)]
    )

    vouchers: Set[str] = set()
    while len(vouchers) < sps_num:
        vouchers.add("".join(random.sample(characters, 3)))

    return sorted(vouchers)



def gene_id_transfer(
    gene2taxa_list: str,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Construct voucher-based gene identifiers and mapping dictionaries.

    Args:
        gene2taxa_list (str): Path to a mapping file with gene IDs in column 1
            and species names in column 2.

    Returns:
        Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, str]]:
            gene2new: Map from original gene IDs to voucher-based gene IDs.
            new2gene: Inverse map from voucher-based gene IDs to originals.
            voucher2taxa: Map from voucher codes to species names.
            taxa2voucher: Map from species names to voucher codes.

    Assumptions:
        The input file contains one gene per row and species names are
        consistent. Voucher assignment is random and therefore non-deterministic
        unless external seeding is applied.
    """
    random.seed(0)
    gene2taxa = read_and_return_dict(gene2taxa_list)
    gene2taxa_sorted = dict(sorted(gene2taxa.items(), key=lambda x: x[1]))

    taxa_list = sorted(set(gene2taxa_sorted.values()))
    taxa2voucher = dict(zip(taxa_list, generate_sps_voucher(len(taxa_list))))
    voucher2taxa = {v: k for k, v in taxa2voucher.items()}

    gene_counter: Dict[str, int] = {}
    new_gene_names: List[str] = []

    for species in gene2taxa_sorted.values():
        gene_counter[species] = gene_counter.get(species, 0) + 1
        new_gene_names.append(f"{taxa2voucher[species]}_{gene_counter[species]}")

    gene2new = dict(zip(gene2taxa_sorted.keys(), new_gene_names))
    new2gene = {v: k for k, v in gene2new.items()}

    return gene2new, new2gene, voucher2taxa, taxa2voucher



def rename_input_tre(tree: Tree, gene2new: Dict[str, str]) -> Tree:
    """Create a renamed copy of a tree using a gene ID mapping.

    Args:
        tree (Tree): ete3 Tree with leaf names as original gene IDs.
        gene2new (Dict[str, str]): Mapping from original gene IDs to new IDs.

    Returns:
        Tree: Copy of the input tree with renamed leaf nodes.

    Assumptions:
        Only leaf names present in gene2new are replaced; others remain unchanged.
    """
    tree_copy = tree.copy()
    for node in tree_copy.traverse():
        if node.name in gene2new:
            node.name = gene2new[node.name]
    return tree_copy


# =========================
# Tree Reading & Rooting
# =========================


def read_tree(tree_path: str) -> Tree:
    """Read a Newick tree as an ete3.Tree, with format fallback.

    Args:
        tree_path (str): Path to a Newick file or a Newick string.

    Returns:
        Tree: Parsed ete3 Tree.

    Assumptions:
        The Newick string is valid; format fallback handles common variants.
    """
    try:
        return Tree(tree_path, format=0)
    except Exception:
        return Tree(tree_path, format=1)



def read_phylo_tree(tree_path: str) -> PhyloTree:
    """Read a Newick tree as an ete3.PhyloTree, with format fallback.

    Args:
        tree_path (str): Path to a Newick file or a Newick string.

    Returns:
        PhyloTree: Parsed ete3 PhyloTree.

    Raises:
        Exception: Propagates parsing errors after logging.

    Assumptions:
        The Newick string is valid; format fallback handles common variants.
    """
    try:
        return PhyloTree(tree_path, format=0)
    except Exception:
        try:
            return PhyloTree(tree_path, format=1)
        except Exception as exc:
            logging.error(f"Failed to read tree: {tree_path}")
            raise exc



def is_rooted(tree: Tree) -> bool:
    """Determine whether a tree is treated as rooted.

    Args:
        tree (Tree): ete3 Tree to evaluate.

    Returns:
        bool: True if the root has exactly two children.

    Assumptions:
        Binary root structure is used as the criterion for rootedness.
    """
    return len(tree.get_children()) == 2



def root_tre_with_midpoint_outgroup(tree: Tree) -> Tree:
    """Root an unrooted tree using midpoint rooting with a safe fallback.

    Args:
        tree (Tree): ete3 Tree to be rooted.

    Returns:
        Tree: A copy of the tree that is rooted if possible.

    Assumptions:
        Trees with fewer than three leaves are returned unchanged because
        midpoint rooting is undefined.
    """
    tree_copy = tree.copy("newick")

    if is_rooted(tree_copy):
        return tree_copy

    leaves = tree_copy.get_leaves()
    if len(leaves) < 3:
        return tree_copy

    try:
        mid = tree_copy.get_midpoint_outgroup()
        tree_copy.set_outgroup(mid if not mid.is_root() else leaves[0])
    except Exception as exc:
        logging.warning(f"Midpoint rooting failed ({exc}); using first leaf")
        tree_copy.set_outgroup(leaves[0])

    return tree_copy


# =========================
# Tree Structure Utilities
# =========================


def num_sptree(sptree):
    idx = 0
    num=1
    for i in sptree.traverse("postorder"):
        i.add_feature("num", num)
        num += 1
        if not i.is_leaf():
            i.name = f"S{idx}"
            idx += 1

    return sptree


def num_tre_node(tree: Tree) -> Tree:
    """Assign deterministic internal-node identifiers in postorder.

    Args:
        tree (Tree): ete3 Tree to label.

    Returns:
        Tree: The same tree with internal nodes named N0, N1, ...

    Assumptions:
        Leaf names are preserved; only internal node names are overwritten.
    """
    idx = 0
    num=1
    for node in tree.traverse("postorder"):
        node.add_feature("num", num)
        num+=1
        if not node.is_leaf():
            node.name = f"N{idx}"
            idx += 1
    return tree



def get_max_deepth(root) -> int:
    """Compute the maximum depth below a node.

    Args:
        root: Tree node with a ``children`` attribute.

    Returns:
        int: Maximum depth from the node to any descendant leaf.

    Assumptions:
        The tree is finite and acyclic.
    """
    if root is None:
        return 0
    return 1 + max((get_max_deepth(c) for c in root.children), default=0)



def calculate_depth(node_a, node_b) -> int:
    """Compute a legacy topological distance between two nodes.

    Args:
        node_a: ete3 node.
        node_b: ete3 node.

    Returns:
        int: Topological distance measured in edges, with a +2 offset when
        one node is an ancestor of the other (legacy behavior).

    Assumptions:
        Nodes belong to the same tree.
    """
    if node_a in node_b.iter_ancestors() or node_b in node_a.iter_ancestors():
        return node_a.get_distance(node_b, topology_only=True) + 2

    ca = node_a.get_common_ancestor(node_b)
    return abs(
        node_a.get_distance(ca, topology_only=True)
        - node_b.get_distance(ca, topology_only=True)
    )



def calculate_deepvar(node_a, node_b) -> int:
    """Compute topological distance between two nodes without legacy offset.

    Args:
        node_a: ete3 node.
        node_b: ete3 node.

    Returns:
        int: Absolute difference in topology-only distances via their LCA.

    Assumptions:
        Nodes belong to the same tree.
    """
    if node_a in node_b.iter_ancestors() or node_b in node_a.iter_ancestors():
        return node_a.get_distance(node_b, topology_only=True)

    ca = node_a.get_common_ancestor(node_b)
    return abs(
        node_a.get_distance(ca, topology_only=True)
        - node_b.get_distance(ca, topology_only=True)
    )


# =========================
# Species Utilities
# =========================


def get_species_list(node) -> List[str]:
    """Extract species labels from leaf names under a node.

    Args:
        node: ete3 node or tree.

    Returns:
        List[str]: Species identifiers parsed as the prefix before the first
        underscore in each leaf name.

    Assumptions:
        Leaf names follow the ``Species_Gene`` convention.
    """
    if node is None:
        return []
    return [leaf.name.split("_")[0] for leaf in node.iter_leaves()]



def get_species_set(node) -> Set[str]:
    """Return unique species labels under a node.

    Args:
        node: ete3 node or tree.

    Returns:
        Set[str]: Unique species identifiers derived from leaf names.

    Assumptions:
        Leaf names follow the ``Species_Gene`` convention.
    """
    return set(get_species_list(node))



def calculate_species_num(node) -> int:
    """Count unique species represented under a node.

    Args:
        node: ete3 node or tree.

    Returns:
        int: Number of unique species under the node.
    """
    return len(get_species_set(node))



def annotate_gene_tree(gene_tree, species_tree):
    """Annotate a gene tree with species mapping and duplication features.

    This function attaches four features to each node: ``map`` (species-tree
    LCA name), ``depth`` (topological depth of the mapped node in the species
    tree), ``overlap`` (species overlap between child clades), and ``is_gd``
    (duplication flag).

    Args:
        gene_tree: ete3 PhyloTree representing a gene family.
        species_tree: ete3 PhyloTree representing the species phylogeny.

    Returns:
        The same ``gene_tree`` with added node features.

    Assumptions:
        Leaf names encode species as the prefix before ``_`` and the species
        tree contains matching leaf labels. Duplication is flagged when the
        overlap between child species sets is at least one (legacy criterion).
        The ``depth`` feature is defined on the species tree (root = 0) and
        reflects the evolutionary depth of the mapped LCA.
    """
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
                overlap_sps = sps_a & sps_b
                overlap_num = len(overlap_sps)

                node.add_feature("overlap", overlap_num)
                node.add_feature("is_gd", overlap_num >= 1)
            else:
                node.add_feature("overlap", 0)
                node.add_feature("is_gd", False)
        else:
            node.add_feature("overlap", 0)
            node.add_feature("is_gd", False)

    return gene_tree


# =========================
# Branch Length Utilities
# =========================


def compute_tip_to_root_branch_length_variance(tree: Tree) -> float:
    """Compute variance of tip-to-root distances.

    Args:
        tree (Tree): ete3 Tree with branch lengths.

    Returns:
        float: Variance of distances from each leaf to the root, or 0.0 for
        fewer than two leaves.

    Assumptions:
        Branch lengths are defined and non-negative.
    """
    distances = [tree.get_distance(leaf) for leaf in tree.iter_leaves()]
    return float(np.var(distances)) if len(distances) > 1 else 0.0



def realign_branch_length(tree: Tree) -> Tree:
    """Reassign branch lengths for visualization-oriented layouts.

    The tree is ladderized and polytomies are resolved; branch lengths are then
    adjusted to align depths for clearer visualization.

    Args:
        tree (Tree): ete3 Tree to modify.

    Returns:
        Tree: The same tree with modified branch lengths.

    Assumptions:
        The tree can be treated as binary after polytomy resolution.
    """
    tree.ladderize()
    tree.resolve_polytomy(recursive=True)
    tree.sort_descendants("support")

    max_depth = get_max_deepth(tree)

    for node in tree.traverse():
        if not node.is_root():
            node.dist = max_depth - get_max_deepth(node) - (
                node.get_distance(tree.get_tree_root()) + 1
            )

    clade_a, clade_b = tree.get_children()
    diff = abs(get_max_deepth(clade_a) - get_max_deepth(clade_b)) + 1
    clade_a.dist += diff
    clade_b.dist += diff

    return tree



def rejust_root_dist(tree: Tree) -> Tree:
    """Adjust root branch distances to balance the two primary clades.

    Args:
        tree (Tree): ete3 Tree whose root branches will be adjusted.

    Returns:
        Tree: The same tree with updated root branch lengths.

    Assumptions:
        The root has exactly two children.
    """
    clade_a, clade_b = tree.get_children()

    def adjust(main, other):
        """Assign branch lengths to emphasize depth balance between clades."""
        main.dist = 1
        other.dist = (
            get_max_deepth(tree) - 1
            if other.is_leaf()
            else get_max_deepth(tree) - get_max_deepth(other)
        )

    if len(clade_a) > len(clade_b):
        adjust(clade_a, clade_b)
    else:
        adjust(clade_b, clade_a)

    return tree


# =========================
# Duplication Detection
# =========================


def judge_support(support: float, threshold: float) -> bool:
    """Evaluate whether a support value meets a threshold.

    Args:
        support (float): Node support, either in [0, 1] or [0, 100].
        threshold (float): Support cutoff, either in [0, 1] or [0, 100].

    Returns:
        bool: True if support >= threshold after normalization.

    Assumptions:
        Values <= 1 are interpreted as fractions and converted to percent.
    """
    if 0 < support < 1:
        support = support * 100
    if 0 < threshold < 1:
        threshold = threshold * 100
    return support >= threshold



def sps_dup_num(sps_list: List[str], unique_sps: Set[str]) -> int:
    """Count the number of duplicated species in a list of species labels.

    Args:
        sps_list (List[str]): Species labels for all leaves in a subtree.
        unique_sps (Set[str]): Unique species set used to initialize counts.

    Returns:
        int: Number of species appearing more than once.

    Assumptions:
        Species labels are consistent and comparable.
    """
    counts = dict.fromkeys(unique_sps, 0)
    duplicated = set()

    for sps in sps_list:
        counts[sps] += 1
        if counts[sps] > 1:
            duplicated.add(sps)

    return len(duplicated)



def get_gene_pairs(gd_node) -> List[Tuple[str, str, str]]:
    """Extract cross-child gene pairs grouped by species for a GD node.

    Args:
        gd_node: Duplication node with exactly two child clades.

    Returns:
        List[Tuple[str, str, str]]: ``(species, gene_left, gene_right)`` tuples.
        Missing genes on one side are represented by ``None``.

    Assumptions:
        Leaf names follow the ``Species_Gene`` convention and each child clade
        corresponds to one duplicated lineage.
    """
    children = gd_node.get_children()
    if len(children) != 2:
        return []

    def collect_genes_by_species(clade):
        species_to_genes = {}
        for leaf in clade.get_leaves():
            species = leaf.name.split("_", 1)[0]
            species_to_genes.setdefault(species, []).append(leaf.name)
        return species_to_genes

    left_map = collect_genes_by_species(children[0])
    right_map = collect_genes_by_species(children[1])
    all_species = set(left_map.keys()) | set(right_map.keys())

    pairs: List[Tuple[str, str, str]] = []
    for species in all_species:
        left_genes = left_map.get(species, [None])
        right_genes = right_map.get(species, [None])
        for gene_left in left_genes:
            for gene_right in right_genes:
                pairs.append((species, gene_left, gene_right))

    return pairs


def find_tre_dup(tree: PhyloTree) -> Tuple[List[str], Set[str]]:
    """Extract duplication event pairs from reconciliation metadata.

    Args:
        tree (PhyloTree): ete3 PhyloTree supporting reconciliation events.

    Returns:
        Tuple[List[str], Set[str]]: A list of duplication event strings and the
        set of leaf names in the tree.

    Raises:
        ValueError: If the tree does not provide reconciliation events.

    Assumptions:
        Events of type ``D`` indicate gene duplications and contain ``in_seqs``
        and ``out_seqs`` attributes.
    """
    if not hasattr(tree, "get_descendant_evol_events"):
        raise ValueError("Tree does not support reconciliation events")

    pairs = []
    leaves = set(tree.get_leaf_names())

    for ev in tree.get_descendant_evol_events():
        if ev.etype == "D":
            pairs.append(",".join(ev.in_seqs) + "<=>" + ",".join(ev.out_seqs))

    return pairs, leaves



def map_species_set_to_node(species_tree, species_set):
    """Map a species set to the corresponding node in a species tree.

    Args:
        species_tree: ete3 PhyloTree for the species phylogeny.
        species_set: Iterable of species names.

    Returns:
        TreeNode or None: The LCA node of the provided species, or None if no
        valid species are found.

    Assumptions:
        Species names correspond to leaf labels in ``species_tree``; missing
        species are ignored.
    """
    if not species_set:
        return None

    nodes = []
    for sp in species_set:
        try:
            nodes.append(species_tree & sp)
        except Exception:
            continue

    if not nodes:
        return None

    # Single species: return the leaf node directly.
    if len(nodes) == 1:
        return nodes[0]

    # Multiple species: return their most recent common ancestor.
    return species_tree.get_common_ancestor(nodes)



def find_dup_node(
    gene_tree: PhyloTree,
    species_tree: PhyloTree,
    gd_support: int = 50,
    clade_support: int = 50,
    dup_species_num: int = 1,
    dup_species_percent: float = 0,
    max_topology_distance: int = 1,
) -> List:
    """Identify duplication nodes using reconciliation and overlap criteria.

    Args:
        gene_tree (PhyloTree): Gene tree annotated with ``map`` and ``is_gd``.
        species_tree (Tree): Species tree used for mapping and distance checks.
        gd_support (int): Minimum support for candidate duplication nodes.
        clade_support (int): Minimum support required for child clades.
        dup_species_num (int): Minimum number of overlapping species required
            to exclude small-scale duplications.
        dup_species_percent (float): Minimum fraction of overlapping species
            relative to all species under the node.
        max_topology_distance (int): Maximum allowed distance between the
            overlap LCA and the mapped duplication node in the species tree.

    Returns:
        List: List of ete3 nodes that meet duplication criteria.

    Assumptions:
        The gene tree has been annotated with ``map`` and ``is_gd`` features and
        uses species labels compatible with ``species_tree``.
    """
    dup_nodes = []

    for node in gene_tree.traverse("postorder"):
        if node.is_leaf() or not getattr(node, "is_gd", False):
            continue
        
        if not judge_support(node.support, gd_support):
            continue

        children = node.get_children()
        if len(children) != 2:
            continue

        if any(not judge_support(c.support, clade_support) for c in children):
            continue

        if not hasattr(node, 'map') or node.map == "unknown":
            continue

        species_set = get_species_set(species_tree&node.map)
        if not species_set:
            continue
        

        # if len(species_set) == 1:
        #     dup_nodes.append(node)
        #     continue

        # if len(species_set) == 2:
        #     sps_tree_lca = species_tree.get_common_ancestor(list(species_set))
        #     s_tree_leaves = set(sps_tree_lca.get_leaf_names())
            
        #     if species_set == s_tree_leaves:
        #         dup_nodes.append(node)
        #         continue

        # else:
        sps_a = get_species_set(children[0])
        sps_b = get_species_set(children[1])
        overlap_sps = sps_a & sps_b
        overlap_num = len(overlap_sps)
        
        # Exclude small-scale duplications by enforcing a minimum overlap.
        if overlap_num < dup_species_num:
            continue

        if overlap_num / len(species_set) < dup_species_percent:
            continue

        # Map overlap species to the species-tree LCA for topology validation.
        dup_map = map_species_set_to_node(species_tree, overlap_sps)

        if dup_map is None:
            continue

        clade_map = species_tree & node.map

        # if abs(children[0].depth-children[1].depth)  >max_topology_distance:
        #     continue

        if abs(clade_map.depth - dup_map.depth) > max_topology_distance:
            continue

        dup_nodes.append(node)

    return dup_nodes



def calculate_gd_num(tree: PhyloTree, species_tree: PhyloTree) -> int:
    """Count gene duplication events using default detection criteria.

    Args:
        tree (PhyloTree): Annotated gene tree.
        species_tree (PhyloTree): Species tree for reconciliation.

    Returns:
        int: Number of duplication nodes.

    Assumptions:
        The tree has been annotated with ``map`` and ``is_gd`` features.
    """
    gd_nodes = find_dup_node(tree, species_tree)
    return len(gd_nodes)
