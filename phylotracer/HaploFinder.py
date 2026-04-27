#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
HaploFinder (refactored) — Haplotype detection and gene conversion analysis.

Key changes vs original
-----------------------
1. Subgenome assignment  — supports diploid, triploid (AAB/ABB) and arbitrary
   polyploids; no longer hard-codes a binary A/B split; distance threshold is
   now topology-aware rather than len(tags).
2. Chromosome homolog pairing — inferred from collinear-block coverage instead
   of the chr1*ratio==chr2 numeric rule.
3. Gene conversion detection — corrected signal polarity (red→blue→red),
   added binomial test + permutation test, requires min_pairs threshold.
4. Gene-pair sorting — 2-D sort on both chromosomes' coordinates to ensure
   genuine collinear ordering before zone scanning.
5. Duplication nodes — handles arbitrary subtree sizes (removes len==3/4 cap),
   emits a WARNING for skipped trees.
"""

import gc
import importlib
import logging
import os
import random
import re
import sys
import time
from collections import Counter, defaultdict
from math import exp, lgamma, log
from typing import Dict, List, Optional, Tuple

import matplotlib
import matplotlib.patches as mpatches
from ete3 import PhyloTree
from tqdm import tqdm

from phylotracer import (
    annotate_gene_tree,
    find_dup_node,
    gene_id_transfer,
    get_species_set,
    judge_support,
    read_and_return_dict,
    read_phylo_tree,
    rename_input_tre,
    root_tre_with_midpoint_outgroup,
)
from phylotracer.Ortho_Retriever import (
    extract_tree,
    iterator,
    offcut_tre,
    rename_OGs_tre_name,
)

matplotlib.use("Agg")
plt = importlib.import_module("matplotlib.pyplot")
plt.rcdefaults()

logger = logging.getLogger(__name__)
logging.getLogger("fontTools").setLevel(logging.WARNING)


# ======================================================
# Section 1: GFF and Lens Parsing  (unchanged interface)
# ======================================================


def read_gff(fn: str) -> Tuple[list, dict]:
    """Read a simplified GFF-like file.

    Returns
    -------
    (data_list, gene_dict)
        gene_dict maps gene_id → [chr, start, end, strand, ?]
    """
    data, data_dict = [], {}
    with open(fn, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split("\t")
            if len(parts) < 6:
                logger.warning(
                    "Skipping malformed GFF row in %s at line %d: "
                    "expected >=6 columns, got %d",
                    fn, line_no, len(parts),
                )
                continue
            gene_id = parts[1]
            gene_info = [parts[0], parts[2], parts[3], parts[4], parts[5]]
            data_dict[gene_id] = gene_info
            data.append(parts)
    return data, data_dict


def read_lens(fn: str, chrs: Optional[str] = None) -> list:
    """Read chromosome length data with optional filtering."""
    data = []
    chrs_lst = None
    if chrs:
        with open(chrs, "r", encoding="utf-8") as f:
            chrs_lst = {i.strip() for i in f if i.strip()}
    with open(fn, "r", encoding="utf-8") as fp:
        for line_no, row in enumerate(fp, start=1):
            stripped = row.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split("\t")
            if len(parts) < 2:
                logger.warning(
                    "Skipping malformed lens row in %s at line %d", fn, line_no
                )
                continue
            if chrs_lst is None or parts[0] in chrs_lst:
                data.append([parts[0], parts[1]])
    return data


# ======================================================
# Section 2: Dotplot Rendering Helpers  (unchanged)
# ======================================================


def _plot_chromosome_axis(lens, gl, gl2, mark, name, is_vertical=False):
    """Unified chromosome axis plotting for both orientations.

    Parameters
    ----------
    is_vertical : bool
        If True, plot vertical axis (chr2 style), else horizontal (chr1 style)
    """
    total_lens = sum([float(k[1]) for k in lens])
    if total_lens == 0:
        return 0
    step = gl / float(total_lens)

    align = dict(
        family="Times New Roman", style="normal",
        horizontalalignment="center", verticalalignment="center",
    )

    if is_vertical:
        # Vertical axis (chr2)
        gl_start, n, start_x = 0.05, 0, 0.95
        mark_y = 0.96
        for k in lens:
            n += float(k[1])
            mark_new = str(mark) + str(k[0])
            x = gl_start + float(n) * step
            mark_x = x - 0.5 * float(k[1]) * step
            plt.plot([x, x], [start_x, start_x - gl2],
                     linestyle="-", color="black", linewidth=0.5)
            plt.text(mark_x, mark_y, mark_new, color="black", fontsize=12,
                     rotation=0, weight="semibold", **align)
        plt.plot([gl_start, gl_start], [start_x, start_x - gl2],
                 linestyle="-", color="black", linewidth=1)
        plt.text(0.5 * (2 * gl_start + gl), mark_y + 0.02, name,
                 color="black", fontsize=18, rotation=0, weight="semibold", **align)
    else:
        # Horizontal axis (chr1)
        gl_start, n, start_x = 0.95, 0, 0.05
        mark_y = 0.04
        for k in lens:
            n += float(k[1])
            mark_new = str(mark) + str(k[0])
            x = gl_start - float(n) * step
            mark_x = x + 0.5 * float(k[1]) * step
            plt.plot([start_x, start_x + gl2], [x, x],
                     linestyle="-", color="black", linewidth=0.5)
            plt.text(mark_y, mark_x, mark_new, color="black", fontsize=12,
                     rotation=90, weight="semibold", **align)
        plt.plot([start_x, start_x + gl2], [gl_start, gl_start],
                 linestyle="-", color="black", linewidth=1)
        plt.text(mark_y - 0.02, 0.5 * (2 * gl_start - gl), name,
                 color="black", fontsize=18, rotation=90, weight="semibold", **align)

    return step


def plot_chr1(lens, gl, gl2, mark, name):
    return _plot_chromosome_axis(lens, gl, gl2, mark, name, is_vertical=False)


def plot_chr2(lens, gl, gl2, mark, name):
    return _plot_chromosome_axis(lens, gl, gl2, mark, name, is_vertical=True)


def gene_location(gff, lens, step):
    loc_gene = {}
    dict_chr = {}
    cumulative = 0
    for chr_name, length in lens:
        dict_chr[chr_name] = cumulative
        cumulative += float(length)

    for k in gff:
        chr_name = k[0]
        if chr_name in dict_chr:
            loc = (dict_chr[chr_name] + float(k[5])) * step
            loc_gene[k[1]] = loc
    return loc_gene


def read_gd_pairs(gd_lst: list) -> dict:
    dict_gd = {}
    for line in gd_lst:
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        gene1, gene2, color = parts[0], parts[1], parts[2]
        pair = f"{gene1}:{gene2}"
        rev_pair = f"{gene2}:{gene1}"
        dict_gd[pair] = color
        if rev_pair not in dict_gd:
            dict_gd[rev_pair] = color
    return dict_gd


def DrawCircle(ax, loc, radius, color, alpha):
    circle = mpatches.Circle(
        loc, radius, edgecolor="none", facecolor=color, alpha=alpha
    )
    ax.add_patch(circle)


def plot_dot(root, loc1, loc2, dict_gd, size=0.001):
    gl_start1, gl_start2 = 0.95, 0.05
    color_order = ["NONE", "red", "green", "purple", "blue"]
    color_map = {
        "NONE": ("gray", 0.6),
        "red": ("red", 0.6),
        "green": ("green", 0.6),
        "purple": ("purple", 0.6),
        "blue": ("blue", 0.6),
    }

    # Pre-group pairs by color to avoid O(5×N) nested loop
    pairs_by_color = defaultdict(list)
    for pair, color in dict_gd.items():
        gene1, gene2 = pair.split(":")
        if gene1 in loc1 and gene2 in loc2:
            pairs_by_color[color].append((gene1, gene2))

    # Draw in color order
    for pass_color in color_order:
        if pass_color not in pairs_by_color:
            continue
        c, a = color_map[pass_color]
        for gene1, gene2 in pairs_by_color[pass_color]:
            x = gl_start1 - loc1[gene1]
            y = gl_start2 + loc2[gene2]
            DrawCircle(root, [y, x], size, c, a)


# ======================================================
# Section 3: Dotplot Pipeline
# ======================================================


def generate_dotplot(
    gff1, gff2, lens1, lens2, gd_pairs, spe1, spe2,
    file_name, target_chr1=None, target_chr2=None, size=None,
    min_pairs: int = 10,
    n_permutations: int = 1000,
    p_threshold: float = 0.05,
):
    """Generate dotplot and run gene-conversion detection."""
    plt.figure(figsize=(10, 10))
    root = plt.axes([0, 0, 1, 1])
    t1 = time.time()
    logger.info("Dotplot of %s starting", file_name)

    gff_1, dict_gff1 = read_gff(gff1)
    gff_2, dict_gff2 = read_gff(gff2)

    if target_chr1 and target_chr2:
        lens_1 = read_lens(lens1, target_chr1)
        lens_2 = read_lens(lens2, target_chr2)
    else:
        lens_1 = read_lens(lens1)
        lens_2 = read_lens(lens2)

    gl1, gl2 = 0.92, 0.92
    step_1 = plot_chr1(lens_1, gl1, gl2, "", spe1)
    step_2 = plot_chr2(lens_2, gl2, gl1, "", spe2)

    dict_gd = read_gd_pairs(gd_pairs)

    # Infer homolog pairs from collinear-block coverage (replaces numeric ratio)
    homolog_pairs = infer_homolog_pairs_from_collinearity(
        dict_gd, dict_gff1, dict_gff2, lens_1, lens_2
    )
    logger.info("Inferred %d homolog chromosome pairs", len(homolog_pairs))

    gene_conversion_list = find_gene_conversion(
        dict_gd, dict_gff1, dict_gff2, homolog_pairs,
        min_pairs=min_pairs,
        n_permutations=n_permutations,
        p_threshold=p_threshold,
    )

    result_conversion = find_gene_pair_info(
        gene_conversion_list, dict_gd, dict_gff1, dict_gff2, file_name
    )
    find_conversion_zones_with_ids_to_file(result_conversion, dict_gff1, dict_gff2)

    # Cache gene locations to avoid redundant GFF traversals
    gene_loc_1 = gene_location(gff_1, lens_1, step_1)
    gene_loc_2 = gene_location(gff_2, lens_2, step_2)

    gc.collect()
    if size:
        plot_dot(root, gene_loc_1, gene_loc_2, dict_gd, size)
    else:
        plot_dot(root, gene_loc_1, gene_loc_2, dict_gd)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    plt.savefig(file_name + "_dotplot.pdf", dpi=500)
    plt.savefig(file_name + "_dotplot.png", dpi=500)
    plt.close()  # Free memory after saving
    logger.info("%s dotplot finished in %.1fs", file_name, time.time() - t1)


# ======================================================
# Section 4: Alignment Labeling  (unchanged interface)
# ======================================================


def assign_colors_by_alignment(
    alignments: Dict[str, List[Tuple[str, str]]],
    alignment_scores: Dict[Tuple[str, str], float]
) -> List[str]:
    """Assign colors to gene pairs based on alignment quality.

    Parameters
    ----------
    alignments : dict
        Mapping from alignment ID to list of (gene_a, gene_b) pairs
    alignment_scores : dict
        Mapping from (gene_a, gene_b) to alignment score

    Returns
    -------
    list
        Sorted list of "gene_a\tgene_b\tcolor" strings
    """
    gene_to_alignments: Dict[str, list] = defaultdict(list)
    for alignment, gene_pairs in alignments.items():
        for gene_a, gene_b in gene_pairs:
            gene_to_alignments[gene_a].append((alignment, gene_b))

    result = set()
    gene_pair_set = set()
    for gene_a, aln_list in gene_to_alignments.items():
        aln_list.sort(key=lambda x: alignment_scores[x[0]], reverse=True)
        for i, (alignment, gene_b) in enumerate(aln_list):
            color = "red" if i == 0 else "blue"
            gene_pair = f"{gene_a}\t{gene_b}"
            if gene_pair not in gene_pair_set:
                result.add(f"{gene_pair}\t{color}")
                gene_pair_set.add(gene_pair)
    return sorted(result)


# ======================================================
# Section 5: Duplication and Orthology Detection
# ======================================================


def generate_combinations(list1: List[str], list2: List[str]) -> List[Tuple[str, str]]:
    return [(e1, e2) for e1 in list1 for e2 in list2]


def get_ortholog_pairs_by_species(t, sp1: str, sp2: str) -> List[Tuple[str, str]]:
    def get_sps_info(node) -> Dict[str, List[str]]:
        sp_count: Dict[str, list] = defaultdict(list)
        for leaf in node.get_leaf_names():
            sp_count[leaf.split("_")[0]].append(leaf)
        return sp_count

    child1, child2 = t.get_children()
    dic1 = get_sps_info(child1)
    dic2 = get_sps_info(child2)

    pairs = (
        generate_combinations(dic1.get(sp1, []), dic2.get(sp2, []))
        + generate_combinations(dic2.get(sp1, []), dic1.get(sp2, []))
    )
    return pairs


def collect_speciation_pairs(
    Phylo_t,
    sp1: str,
    sp2: str,
    new_named_gene2gene_dic: dict,
    processed_lines: list,
    written_results: set,
    pair_support: int = 50,
) -> None:
    """Collect speciation-based ortholog pairs from phylogenetic tree.

    Parameters
    ----------
    Phylo_t : PhyloTree
        Annotated phylogenetic tree
    sp1, sp2 : str
        Species identifiers to compare
    new_named_gene2gene_dic : dict
        Mapping from renamed genes to original gene IDs
    processed_lines : list
        Output list to append gene pairs
    written_results : set
        Set to track already written pairs
    pair_support : int
        Minimum support value for pair nodes
    """
    events = Phylo_t.get_descendant_evol_events()
    for event in events:
        if event.etype != "S":
            continue
        event_nodes = set(event.in_seqs) | set(event.out_seqs)
        spec_node = Phylo_t.get_common_ancestor(list(event_nodes))
        sp_set = get_species_set(spec_node)
        if len(sp_set) != 2 or sp1 not in sp_set or sp2 not in sp_set:
            continue
        for gene1, gene2 in get_ortholog_pairs_by_species(spec_node, sp1, sp2):
            pair_node = Phylo_t.get_common_ancestor([gene1, gene2])
            if judge_support(pair_node.support, pair_support):
                gene_a = new_named_gene2gene_dic[gene1]
                gene_b = new_named_gene2gene_dic[gene2]
                token = f"{gene2}-r"
                if token not in written_results:
                    processed_lines.append(f"{gene_a}\t{gene_b}\tred\n")
                    written_results.add(token)


def collect_duplication_nodes_legacy(
    Phylo_t,
    sp1: str,
    sp2: str,
    support_value: int = 50
) -> list:
    """Collect duplication nodes from phylogenetic tree.

    Parameters
    ----------
    Phylo_t : PhyloTree
        Annotated phylogenetic tree
    sp1, sp2 : str
        Species identifiers to filter
    support_value : int
        Minimum support threshold for nodes

    Returns
    -------
    list
        List of duplication nodes meeting criteria
    """
    dup_node_list = []
    events = Phylo_t.get_descendant_evol_events()
    for event in events:
        if event.etype != "D":
            continue
        event_nodes = set(event.in_seqs) | set(event.out_seqs)
        dup_node = Phylo_t.get_common_ancestor(list(event_nodes))
        if not judge_support(dup_node.support, support_value):
            continue
        children = dup_node.get_children()
        if len(children) != 2:
            continue
        child1, child2 = children
        if not (
            judge_support(child1.support, support_value)
            and judge_support(child2.support, support_value)
        ):
            continue
        sp_set = get_species_set(dup_node)
        if sp1 in sp_set and sp2 in sp_set:
            dup_node_list.append(dup_node)
    return dup_node_list


def find_independent_dup_nodes(dup_node_list: list) -> list:
    dup_node_list.sort(key=lambda x: len(x.get_leaf_names()))
    node_dic: dict = {}
    for node in dup_node_list:
        if len(get_species_set(node)) == 1:
            continue
        tips = len(node)
        node_dic.setdefault(tips, []).append(node)
    if node_dic:
        return min(node_dic.items(), key=lambda x: x[0])[1]
    return []


# ======================================================
# Section 6a: Chromosome Homolog Pair Inference
#             Replaces the hard-coded chr1*ratio==chr2 rule
# ======================================================


def infer_homolog_pairs_from_collinearity(
    dict_gd: dict,
    dict_gff1: dict,
    dict_gff2: dict,
    lens_1: list,
    lens_2: list,
    min_shared_pairs: int = 5,
) -> List[Tuple[str, str]]:
    """Infer chromosome homolog pairs from collinear-block gene-pair coverage.

    Method
    ------
    For every (chr_a, chr_b) combination, count how many gene pairs from
    dict_gd map to that combination.  For each chr_a, the chr_b that receives
    the most gene pairs is its primary homolog (and vice-versa).  Pairs that
    do not meet min_shared_pairs on both directions are discarded.

    This is a coverage-based proxy for synteny; it requires no assumption
    about chromosome-naming conventions.

    Parameters
    ----------
    dict_gd : dict
        Gene-pair → color mapping (output of read_gd_pairs).
    dict_gff1, dict_gff2 : dict
        Gene → [chr, start, end, …] mappings from read_gff.
    lens_1, lens_2 : list
        [[chr, length], …] lists from read_lens.
    min_shared_pairs : int
        Minimum gene pairs required for a pair to be considered homologous.

    Returns
    -------
    List of (chr_a, chr_b) tuples representing inferred homolog pairs.
    """
    chrs_1 = {row[0] for row in lens_1}
    chrs_2 = {row[0] for row in lens_2}

    # Count gene pairs per (chr_a, chr_b) combination
    pair_count: Dict[Tuple[str, str], int] = defaultdict(int)
    for pair_key in dict_gd:
        gene_a, gene_b = pair_key.split(":", 1)
        gff_a = dict_gff1.get(gene_a)
        gff_b = dict_gff2.get(gene_b)
        if gff_a and gff_b:
            chr_a, chr_b = gff_a[0], gff_b[0]
            if chr_a in chrs_1 and chr_b in chrs_2:
                pair_count[(chr_a, chr_b)] += 1

    # For each chr_a: best chr_b by count
    best_b_for_a: Dict[str, Tuple[str, int]] = {}
    for (chr_a, chr_b), cnt in pair_count.items():
        if chr_a not in best_b_for_a or cnt > best_b_for_a[chr_a][1]:
            best_b_for_a[chr_a] = (chr_b, cnt)

    # For each chr_b: best chr_a by count
    best_a_for_b: Dict[str, Tuple[str, int]] = {}
    for (chr_a, chr_b), cnt in pair_count.items():
        if chr_b not in best_a_for_b or cnt > best_a_for_b[chr_b][1]:
            best_a_for_b[chr_b] = (chr_a, cnt)

    # Reciprocal best hits + minimum coverage filter
    homolog_pairs = []
    for chr_a, (chr_b, cnt_ab) in best_b_for_a.items():
        if cnt_ab < min_shared_pairs:
            continue
        if best_a_for_b.get(chr_b, (None,))[0] == chr_a:
            homolog_pairs.append((chr_a, chr_b))
            logger.info(
                "Homolog pair inferred: %s <-> %s  (shared pairs=%d)",
                chr_a, chr_b, cnt_ab,
            )

    if not homolog_pairs:
        logger.warning(
            "No homolog chromosome pairs could be inferred "
            "(min_shared_pairs=%d). Consider lowering --min_shared_pairs.",
            min_shared_pairs,
        )
    return homolog_pairs


# ======================================================
# Section 6b: Gene Conversion Detection  (corrected)
# ======================================================


def _binomial_p_value(blue_count: int, total: int, background_blue_rate: float) -> float:
    """One-tailed binomial P(X >= blue_count | n=total, p=background_blue_rate).

    Uses a log-space PMF seed plus recurrence to avoid float overflow for
    large ``total``.
    """
    if total == 0 or background_blue_rate <= 0:
        return 1.0
    if blue_count <= 0:
        return 1.0
    if blue_count > total:
        return 0.0
    if background_blue_rate >= 1:
        return 1.0 if blue_count <= total else 0.0

    p = background_blue_rate
    k = blue_count

    log_pmf_k = (
        lgamma(total + 1)
        - lgamma(k + 1)
        - lgamma(total - k + 1)
        + k * log(p)
        + (total - k) * log(1 - p)
    )
    pmf = exp(log_pmf_k)
    tail = pmf

    for j in range(k, total):
        pmf *= ((total - j) / (j + 1)) * (p / (1 - p))
        tail += pmf

    return max(0.0, min(1.0, tail))


def _permutation_p_value(
    blue_count: int,
    total: int,
    all_colors: list,
    n_permutations: int = 1000,
    rng_seed: int = 42,
) -> float:
    """Permutation test: how often does random sampling yield >= blue_count blues?"""
    if total == 0:
        return 1.0
    rng = random.Random(rng_seed)
    exceed = 0
    for _ in range(n_permutations):
        sample = rng.choices(all_colors, k=total)
        if sample.count("blue") >= blue_count:
            exceed += 1
    return exceed / n_permutations


def find_gene_conversion(
    dict_gd: dict,
    dict_gff1: dict,
    dict_gff2: dict,
    homolog_pairs: List[Tuple[str, str]],
    min_pairs: int = 10,
    n_permutations: int = 1000,
    p_threshold: float = 0.05,
) -> list:
    """Collect homolog chromosome pairs to scan for local conversion islands.

    Gene conversion is defined downstream as a local continuous blue interval
    interrupting a red collinear background (red→blue→red). Therefore this
    function deliberately does not require chromosome-wide blue enrichment:
    such a filter would discard the expected red-dominant chromosome pairs
    before local islands can be scanned.

    Parameters
    ----------
    homolog_pairs : list of (chr_a, chr_b)
        Output of infer_homolog_pairs_from_collinearity.
    min_pairs : int
        Minimum number of gene pairs on a chromosome pair to scan.
    n_permutations : int
        Retained for CLI compatibility; not used by the local-island scanner.
    p_threshold : float
        Retained for CLI compatibility; not used by the local-island scanner.
    """
    all_colors = [c for c in dict_gd.values() if c in ("red", "blue")]
    total_global = len(all_colors)
    bg_blue_rate = all_colors.count("blue") / total_global if total_global > 0 else 0.0
    logger.info(
        "Genome-wide blue rate: %.3f  (total pairs=%d)", bg_blue_rate, total_global
    )

    # Pre-build chromosome lookup for faster filtering
    chr_pairs_lookup: Dict[Tuple[str, str], List[Tuple[str, str, str]]] = defaultdict(list)
    for pair_key, color in dict_gd.items():
        if color not in ("red", "blue"):
            continue
        gene_a, gene_b = pair_key.split(":", 1)
        gff_a = dict_gff1.get(gene_a)
        gff_b = dict_gff2.get(gene_b)
        if gff_a and gff_b:
            chr_a, chr_b = gff_a[0], gff_b[0]
            chr_pairs_lookup[(chr_a, chr_b)].append((gene_a, gene_b, color))

    block_list = []
    for chr_a, chr_b in homolog_pairs:
        pairs = chr_pairs_lookup.get((chr_a, chr_b), [])
        red_count = sum(1 for _, _, c in pairs if c == "red")
        blue_count = sum(1 for _, _, c in pairs if c == "blue")
        total = red_count + blue_count

        if total < min_pairs:
            logger.debug(
                "Skipping %s<->%s: only %d pairs (min_pairs=%d)",
                chr_a, chr_b, total, min_pairs,
            )
            continue

        red_ratio = red_count / total
        blue_ratio = blue_count / total

        block = dict(
            chr_a=chr_a, chr_b=chr_b,
            red_count=red_count, blue_count=blue_count,
            red_ratio=red_ratio, blue_ratio=blue_ratio,
        )
        logger.info(
            "%s<->%s  red=%d blue=%d  blue_ratio=%.3f  queued for red-blue-red scan",
            chr_a, chr_b, red_count, blue_count, blue_ratio,
        )
        block_list.append(block)

    return block_list


def find_gene_pair_info(
    gene_conversion_list: list,
    dict_gd: dict,
    dict_gff1: dict,
    dict_gff2: dict,
    gd_pairs: str,
) -> list:
    """Collect and sort gene pairs for significant conversion blocks.

    Fix: 2-D sort on (chr_a, start_a, chr_b, start_b) instead of single-axis.
    This ensures collinear ordering on both chromosomes before zone scanning.
    """
    # Build set of target chromosome pairs for faster filtering
    target_chr_pairs = {(block["chr_a"], block["chr_b"]) for block in gene_conversion_list}

    sort_lst = []
    for pair_key, color in dict_gd.items():
        gene_a, gene_b = pair_key.split(":", 1)
        info_a = dict_gff1.get(gene_a)
        info_b = dict_gff2.get(gene_b)
        if not info_a or not info_b:
            continue

        chr_a, chr_b = info_a[0], info_b[0]
        if (chr_a, chr_b) not in target_chr_pairs:
            continue

        try:
            start_a = int(info_a[1])
            start_b = int(info_b[1])
        except (TypeError, ValueError):
            continue
        sort_lst.append((gene_a, gene_b, color, chr_a, start_a, chr_b, start_b))

    # 2-D sort: primary on chr_a position, secondary on chr_b position
    sort_lst.sort(key=lambda x: (x[3], x[4], x[5], x[6]))

    new_lst = []
    with open(f"gene_conversion_{gd_pairs}.txt", "w") as f:
        for gene_a, gene_b, color, *_ in sort_lst:
            f.write(f"{gene_a}\t{gene_b}\t{color}\n")
            new_lst.append((gene_a, gene_b, color))
    return new_lst


def find_conversion_zones_with_ids_to_file(
    data: list,
    dict_gff1: dict,
    dict_gff2: dict,
    output_file: str = "gene_conversion.txt",
) -> None:
    """Scan for gene conversion zones using the corrected red→blue→red pattern.

    Original code scanned blue→red→blue, which is the inverted signal.

    Gene conversion causes a donor-sequence island to interrupt the expected
    collinear (red) background, producing a local enrichment of blue (cross-clade)
    pairs flanked by red (in-clade) pairs on both sides.

    Corrected pattern
    -----------------
        red … red   |   blue … blue   |   red … red
        (normal)        (converted)        (normal)
    """

    def _consume_color_run(data, start_idx, expected_color, chrom1, chrom2, dict_gff1, dict_gff2):
        """Consume consecutive pairs of the same color on the same chromosome pair.

        Returns (end_idx, run_length) where end_idx is the last valid index.
        """
        i = start_idx
        n = len(data)
        while i < n:
            g1, g2, c = data[i]
            if c != expected_color:
                break
            gff1 = dict_gff1.get(g1)
            gff2 = dict_gff2.get(g2)
            if not gff1 or not gff2 or gff1[0] != chrom1 or gff2[0] != chrom2:
                break
            i += 1
        return i - 1, i - start_idx

    conversion_zones = []
    n = len(data)
    zone_id = 1
    i = 0

    while i < n - 1:
        gene1_start, gene2_start, color = data[i]

        # Must start with a red run
        if color != "red":
            i += 1
            continue

        gff1_info = dict_gff1.get(gene1_start)
        gff2_info = dict_gff2.get(gene2_start)
        if not gff1_info or not gff2_info:
            i += 1
            continue

        chrom1, chrom2 = gff1_info[0], gff2_info[0]

        # Consume leading red run
        start_red1 = i
        end_red1, _ = _consume_color_run(data, i, "red", chrom1, chrom2, dict_gff1, dict_gff2)
        i = end_red1 + 1

        if i >= n or data[i][2] != "blue":
            continue

        # Consume blue (converted) run
        start_blue = i
        end_blue, blue_run_len = _consume_color_run(data, i, "blue", chrom1, chrom2, dict_gff1, dict_gff2)
        i = end_blue + 1

        if i >= n or data[i][2] != "red":
            i = start_blue
            continue

        # Consume trailing red run
        start_red2 = i
        end_red2, _ = _consume_color_run(data, i, "red", chrom1, chrom2, dict_gff1, dict_gff2)
        i = end_red2 + 1

        # Require at least 3 blue pairs in the converted run
        if blue_run_len >= 3:
            # Report the local island only: nearest red anchor on each side
            # plus the continuous blue run. The full red runs are background
            # context and would make candidate zones artificially large.
            local_start = end_red1
            local_end = start_red2
            conversion_zones.append((zone_id, local_start, local_end,
                                     start_blue, end_blue, blue_run_len))
            logger.info(
                "Conversion zone %d: %s&%s  blue_run=%d  [idx %d-%d]",
                zone_id, chrom1, chrom2, blue_run_len, start_blue, end_blue,
            )
            zone_id += 1

        # Allow overlapping detection from start_red2
        i = start_red2

    with open(output_file, "w") as f:
        f.write(
            "# zone_id\tcontig1\tcontig2\t"
            "gene1\tstart1\tend1\tgene2\tstart2\tend2\tcolor\n"
        )
        for zone in conversion_zones:
            zid, s_r1, e_r2, s_b, e_b, blen = zone
            first_gene1, first_gene2, _ = data[s_r1]
            contig1 = dict_gff1[first_gene1][0]
            contig2 = dict_gff2[first_gene2][0]
            f.write(
                f"# Conversion Zone {zid}: {contig1}&{contig2} "
                f"(blue_run={blen} pairs)\n"
            )
            for j in range(s_r1, e_r2 + 1):
                if j >= n:
                    break
                g1, g2, color = data[j]
                gff1 = dict_gff1.get(g1)
                gff2 = dict_gff2.get(g2)
                if gff1 and gff2:
                    c1, ps1, pe1 = gff1[0], gff1[1], gff1[2]
                    c2, ps2, pe2 = gff2[0], gff2[1], gff2[2]
                    f.write(
                        f"{zid}\t{c1}\t{c2}\t"
                        f"{g1}\t{ps1}\t{pe1}\t{g2}\t{ps2}\t{pe2}\t{color}\n"
                    )

    logger.info("Wrote %d conversion zones to %s", len(conversion_zones), output_file)


# ======================================================
# Section 7: GD Pair Processing  (duplication node cap removed)
# ======================================================


def process_gd_result(gf, imap, input_sps_tree, sp1, sp2, support, pair_support):
    """Process gene-family trees to assign red/blue labels.

    Fix: removed the len(i)==3 / len(i)==4 cap.  All duplication nodes with
    ≥3 leaves and both target species present are now processed.  Trees that
    are skipped due to missing species are logged at DEBUG level.
    """
    tre_dic = read_and_return_dict(gf)
    (gene2new_named_gene_dic, new_named_gene2gene_dic,
     voucher2taxa_dic, taxa2voucher_dic) = gene_id_transfer(imap)
    rename_sp1 = taxa2voucher_dic[sp1]
    rename_sp2 = taxa2voucher_dic[sp2]
    renamed_sptree = rename_input_tre(
        read_phylo_tree(input_sps_tree), taxa2voucher_dic
    )
    processed_lines: list = []
    written_results: set = set()
    skipped_trees = 0

    for tre_ID, tre_path in tre_dic.items():
        Phylo_t0 = read_phylo_tree(tre_path)
        Phylo_t0.resolve_polytomy(recursive=True)
        Phylo_t0.sort_descendants()
        Phylo_t1 = rename_input_tre(Phylo_t0, gene2new_named_gene_dic)

        sps_tol = get_species_set(Phylo_t1)
        if len(sps_tol) == 1:
            skipped_trees += 1
            continue
        if len(sps_tol) == 2 and not {rename_sp1, rename_sp2}.issubset(sps_tol):
            skipped_trees += 1
            continue

        annotate_gene_tree(Phylo_t1, renamed_sptree)
        dup_node_list = collect_duplication_nodes_legacy(
            Phylo_t1, rename_sp1, rename_sp2, support_value=support
        )
        collect_speciation_pairs(
            Phylo_t1, rename_sp1, rename_sp2,
            new_named_gene2gene_dic, processed_lines, written_results, pair_support,
        )

        for dup_node in dup_node_list:
            n_leaves = len(dup_node.get_leaf_names())
            if n_leaves < 3:
                continue  # degenerate node; skip silently

            gd_clade1, gd_clade2 = dup_node.get_children()
            tips1 = set(gd_clade1.get_leaf_names())
            tips2 = set(gd_clade2.get_leaf_names())
            sp_set1 = get_species_set(gd_clade1)
            sp_set2 = get_species_set(gd_clade2)

            # Both clades must contain at least one of the two target species
            if not (
                (rename_sp1 in sp_set1 or rename_sp2 in sp_set1)
                and (rename_sp1 in sp_set2 or rename_sp2 in sp_set2)
            ):
                logger.debug(
                    "Tree %s dup node (n=%d): species not split across clades, skipping",
                    tre_ID, n_leaves,
                )
                continue

            gene_pairs = get_ortholog_pairs_by_species(
                dup_node, rename_sp1, rename_sp2
            )
            for gene1, gene2 in gene_pairs:
                item_set = {gene1, gene2}
                gene_a = new_named_gene2gene_dic[gene1]
                gene_b = new_named_gene2gene_dic[gene2]
                if item_set <= tips1 or item_set <= tips2:
                    color, suffix = "red", "r"
                else:
                    color, suffix = "blue", "b"
                token = f"{gene2}-{suffix}"
                if token not in written_results:
                    processed_lines.append(f"{gene_a}\t{gene_b}\t{color}\n")
                    written_results.add(token)

    if skipped_trees:
        logger.info("Skipped %d trees (single-species or missing target sp)", skipped_trees)

    sorted_lines = sorted(processed_lines, key=lambda x: x.split("\t")[0])
    with open("color_label.txt", "w") as fh:
        for line in sorted_lines:
            fh.write(line)
    return sorted_lines


# ======================================================
# Section 8: FASTA Utilities  (unchanged)
# ======================================================


def read_fasta(file_path: str) -> dict:
    sequences: dict = {}
    current_id = None
    current_seq: list = []
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences


def write_fasta(sequences: dict, file_path: str) -> None:
    with open(file_path, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")


# ======================================================
# Section 9: Single-Copy Tree Extraction  (unchanged)
# ======================================================


def get_single_copy_tree(tree, tree_id, tree_path, gene2new_named_gene_dic,
                         renamed_length_dict=None):
    renamed_len_dic = {k: 0 for k in tree.get_leaf_names()}
    principal_gene_set, filtered_offcut_ev_seqs = offcut_tre(tree, renamed_len_dic)
    minor_orthologs: list = []
    minor_orthologs = iterator(
        filtered_offcut_ev_seqs, tree, gene2new_named_gene_dic,
        minor_orthologs, tree_path, renamed_length_dict,
    )
    ordered_name_OG_list = rename_OGs_tre_name(
        principal_gene_set, minor_orthologs, tree_id
    )
    sctree = []
    for tre_name, OG_set in ordered_name_OG_list:
        phylo_tree_0 = read_phylo_tree(tree_path)
        phylo_tree = root_tre_with_midpoint_outgroup(phylo_tree_0)
        rename_t = rename_input_tre(phylo_tree, gene2new_named_gene_dic)
        sctree.append(extract_tree(OG_set, rename_t))
    return sctree


# ======================================================
# Section 10: Hybrid Subgenome Assignment  (polyploid-aware)
# ======================================================


def assign_hybrid_subgenome(
    tree,
    hybrid_prefix: str,
    diploid_tags: List[str],
) -> Dict[str, str]:
    """Assign each hybrid gene to a parental subgenome.

    Each parental species tag maps to one stable label (A, B, C, ...). For each
    hybrid gene, assignment is based on the closest parental gene in the gene
    tree topology rather than chromosome naming or tree-local copy numbering.

    Parameters
    ----------
    tree : PhyloTree
        Gene tree (already renamed).
    hybrid_prefix : str
        String present in all hybrid gene names.
    diploid_tags : list of str
        Parental species tags, e.g. ['Ath', 'Aly'] for diploid or
        ['Pa', 'Pb', 'Pc'] for triploid.
    """
    subgenome_base_labels = {tag: chr(ord("A") + i)
                              for i, tag in enumerate(diploid_tags)}
    result: Dict[str, str] = {}

    # Pre-build sets for faster lookups
    diploid_tags_set = set(diploid_tags)

    def is_hybrid(name: str) -> bool:
        return hybrid_prefix in name

    def get_diploid_tag(name: str) -> Optional[str]:
        for tag in diploid_tags_set:
            if tag in name:
                return tag
        return None

    def _label_for_tag(tag: str) -> str:
        return subgenome_base_labels[tag]

    def sister_parent_tag(leaf_node) -> Optional[str]:
        """Find the nearest sister clade carrying one unambiguous parent tag."""
        child = leaf_node
        node = leaf_node.up
        while node is not None:
            sister_tags = set()
            sister_has_hybrid = False
            for sister in node.children:
                if sister is child:
                    continue
                for leaf_name in sister.get_leaf_names():
                    if is_hybrid(leaf_name):
                        sister_has_hybrid = True
                    tag = get_diploid_tag(leaf_name)
                    if tag:
                        sister_tags.add(tag)

            if len(sister_tags) == 1 and not sister_has_hybrid:
                return next(iter(sister_tags))
            if len(sister_tags) > 1:
                return None

            child = node
            node = node.up
        return None

    for leaf_node in tree.iter_leaves():
        leaf_name = leaf_node.name
        if not is_hybrid(leaf_name):
            continue
        tag = sister_parent_tag(leaf_node)
        if tag is None:
            result[leaf_name] = "unknown"
            logger.debug(
                "Unresolved hybrid leaf by nearest-sister topology: %s",
                leaf_name,
            )
        else:
            result[leaf_name] = _label_for_tag(tag)

    return result


# ======================================================
# Section 11: Sequence Splitting
# ======================================================


def split_sequences(
    input_GF_list, input_imap, hyb_sps, parental_sps,
    gff, input_fasta, cluster_file, chrs_per_subgenome=10,
):
    """Split hybrid sequences into per-subgenome FASTA files.

    For polyploids with >2 parental species, additional subgenome files
    (C, D, …) are written automatically based on subgenome labels present
    in the assignment output.
    """
    tre_dic = read_and_return_dict(input_GF_list)
    gff_1, dict_gff1 = read_gff(gff)
    _ = gff_1
    (gene2new_named_gene_dic, new_named_gene2gene_dic,
     voucher2taxa_dic, taxa2voucher_dic) = gene_id_transfer(input_imap)
    hyb_sps_renamed = taxa2voucher_dic.get(hyb_sps, hyb_sps)
    parental_sps_renamed = [taxa2voucher_dic.get(s, s) for s in parental_sps]

    all_gene_labels: dict = {}
    conflict_genes: set = set()

    logger.info(
        "split mode: cluster_file=%s (used for logging only)", cluster_file
    )

    for tre_ID, tre_path in tre_dic.items():
        Phylo_t0 = read_phylo_tree(tre_path)
        Phylo_t0.resolve_polytomy(recursive=True)
        Phylo_t0.sort_descendants()
        Phylo_t1 = rename_input_tre(Phylo_t0, gene2new_named_gene_dic)
        leaf_names = Phylo_t1.get_leaf_names()
        leaf_species = {g.split("_")[0] for g in leaf_names}
        if len(leaf_species) == 1:
            continue
        if not any(hyb_sps_renamed in ln for ln in leaf_names):
            continue

        a = assign_hybrid_subgenome(Phylo_t1, hyb_sps_renamed, parental_sps_renamed)
        b = {new_named_gene2gene_dic[k]: v for k, v in a.items()}

        for k1, v1 in b.items():
            if k1 in all_gene_labels and all_gene_labels[k1] != v1:
                all_gene_labels[k1] = "unknown"
                conflict_genes.add(k1)
            else:
                all_gene_labels[k1] = v1

    all_sequences = read_fasta(input_fasta)
    out_dir = "haplofinder_split"
    os.makedirs(out_dir, exist_ok=True)

    # Collect all distinct subgenome labels (supports A, B, C, A_2, etc.)
    all_labels = set(all_gene_labels.values()) - {"unknown"}
    split_buckets: Dict[str, dict] = {lbl: {} for lbl in all_labels}
    split_buckets["unknown"] = {}
    assigned_total = 0
    skipped_unassigned = 0

    with open(os.path.join(out_dir, "split_assignment.tsv"), "w") as out:
        out.write("gene_id\tsubgenome\tstatus\tchromosome\n")
        for gene_id, seq in all_sequences.items():
            label = all_gene_labels.get(gene_id)
            if label is None:
                skipped_unassigned += 1
                continue
            assigned_total += 1
            final_label = "unknown" if gene_id in conflict_genes else label
            status = "conflict" if gene_id in conflict_genes else "ok"
            chr_name = dict_gff1.get(gene_id, [None])[0]
            split_buckets.setdefault(final_label, {})[gene_id] = seq
            out.write(
                f"{gene_id}\t{final_label}\t{status}\t"
                f"{chr_name if chr_name else 'NA'}\n"
            )

    for lbl, seqs in split_buckets.items():
        safe_lbl = lbl.replace("/", "_")
        write_fasta(seqs, os.path.join(out_dir, f"split_subgenome_{safe_lbl}.fasta"))

    with open(os.path.join(out_dir, "split_summary.txt"), "w") as out:
        out.write(f"input_fasta_total={len(all_sequences)}\n")
        out.write(f"assigned_total={assigned_total}\n")
        out.write(f"skipped_unassigned={skipped_unassigned}\n")
        for lbl, seqs in split_buckets.items():
            out.write(f"subgenome_{lbl}={len(seqs)}\n")
        out.write(f"conflict_genes={len(conflict_genes)}\n")

    logger.info("split mode outputs written to %s", os.path.abspath(out_dir))


# ======================================================
# Section 12: Chromosome Mapping (kept for backward compat)
# ======================================================


def get_chromosome_subgenome(chr_name: str, chrs_per_subgenome: int = 10) -> Optional[str]:
    """Legacy numeric mapping; used only as a validation helper now.

    Do NOT use for primary subgenome assignment — see assign_hybrid_subgenome.
    """
    match = re.search(r"(\d+)", chr_name)
    if match:
        n = int(match.group(1))
        if 1 <= n <= chrs_per_subgenome:
            return "A"
        elif chrs_per_subgenome + 1 <= n <= 2 * chrs_per_subgenome:
            return "B"
    return None


# ======================================================
# Section 13: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="HaploFinder (refactored) — polyploid subgenome assignment "
                    "and gene conversion detection"
    )
    parser.add_argument("gff1")
    parser.add_argument("gff2")
    parser.add_argument("lens1")
    parser.add_argument("lens2")
    parser.add_argument("--blastp_pairs", default=None,
                        help="Optional: path to BLAST pair file for additional coloring")
    parser.add_argument("--synteny_pairs", default=None,
                        help="Optional: path to synteny pairs file for alignment-based coloring")
    parser.add_argument("spe1")
    parser.add_argument("spe2")
    parser.add_argument("num", type=int)
    parser.add_argument("gf")
    parser.add_argument("imap")
    parser.add_argument("target_chr1")
    parser.add_argument("target_chr2")
    parser.add_argument("size", type=float)
    parser.add_argument("input_sps_tree")
    parser.add_argument("--min_shared_pairs", type=int, default=5,
                        help="Min gene pairs to infer a chromosome homolog pair")
    parser.add_argument("--min_conv_pairs", type=int, default=10,
                        help="Min gene pairs on a chr pair to test gene conversion")
    parser.add_argument("--n_permutations", type=int, default=1000)
    parser.add_argument("--p_threshold", type=float, default=0.05)
    args = parser.parse_args()

    process_blastp_pairs = _process_blastp_result_stub(args.blastp_pairs, args.num)
    alignments, alignment_scores = _parse_synteny_file_stub(args.synteny_pairs)
    process_synteny_pairs = assign_colors_by_alignment(alignments, alignment_scores)
    process_gd_pairs = process_gd_result(
        args.gf, args.imap, args.input_sps_tree, args.spe1, args.spe2, 50, 50
    )

    for pairs, tag in [
        (process_blastp_pairs, "blastp_pairs"),
        (process_synteny_pairs, "synteny_pairs"),
        (process_gd_pairs, "gd_pairs"),
    ]:
        generate_dotplot(
            args.gff1, args.gff2, args.lens1, args.lens2,
            pairs, args.spe1, args.spe2, tag,
            args.target_chr1, args.target_chr2, args.size,
            min_pairs=args.min_conv_pairs,
            n_permutations=args.n_permutations,
            p_threshold=args.p_threshold,
        )

    total_lst = _process_total_color_list_stub(
        [process_blastp_pairs, process_synteny_pairs, process_gd_pairs]
    )
    generate_dotplot(
        args.gff1, args.gff2, args.lens1, args.lens2,
        total_lst, args.spe1, args.spe2, "total_pairs",
        args.target_chr1, args.target_chr2, args.size,
        min_pairs=args.min_conv_pairs,
        n_permutations=args.n_permutations,
        p_threshold=args.p_threshold,
    )


# ============================================================
# Optional input processors (fallback implementations)
# ============================================================

def _process_blastp_result_stub(blastp_pairs, num):
    """Fallback: read simple TSV with chr1,s1,e1,chr2,s2,e2,color format."""
    if blastp_pairs is None or not os.path.exists(blastp_pairs):
        logger.info("BLAST pairs file not provided, skipping blast-based coloring")
        return []
    result = []
    with open(blastp_pairs, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                try:
                    result.append((parts[0], int(parts[1]), parts[2],
                                    parts[3], int(parts[4]), parts[5]))
                except ValueError:
                    continue
    return result


def _parse_synteny_file_stub(synteny_pairs):
    """Fallback: parse simple TSV for alignment-based coloring."""
    if synteny_pairs is None or not os.path.exists(synteny_pairs):
        logger.info("Synteny pairs file not provided, skipping synteny-based coloring")
        return {}, {}
    alignments = {}
    scores = {}
    with open(synteny_pairs, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                aln_id = parts[0]
                gene_a, gene_b = parts[1].split(",") if "," in parts[1] else (parts[1], "")
                alignments.setdefault(aln_id, []).append((gene_a.strip(), gene_b.strip()))
                try:
                    scores[(gene_a.strip(), gene_b.strip())] = float(parts[2])
                except ValueError:
                    scores[(gene_a.strip(), gene_b.strip())] = 1.0
    return alignments, scores


def _process_total_color_list_stub(total_pairs):
    """Merge multiple pair lists into one (used for total summary plot)."""
    merged = []
    for pair_list in total_pairs:
        merged.extend(pair_list)
    return merged
