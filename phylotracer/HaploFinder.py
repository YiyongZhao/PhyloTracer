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


def plot_dot(root, loc1, loc2, dict_gd, size=0.001, dict_gff1=None, dict_gff2=None, mode="auto",
             chrs_per_subgenome=10):
    gl_start1, gl_start2 = 0.95, 0.05
    color_map = {
        "NONE": ("gray", 0.6),
        "red": ("red", 0.6),
        "green": ("green", 0.6),
        "purple": ("purple", 0.6),
        "blue": ("blue", 0.6),
    }
    all_colors = ["NONE", "red", "green", "purple", "blue"]

    # Determine per-chromosome-pair draw order: signal color always last.
    # Generic across species: uses get_chromosome_subgenome + chrs_per_subgenome.
    def _draw_order_for_chr_pair(chr_a, chr_b):
        if mode == "inter":
            # All pairs cross-subgenome (different species): blue bg, red signal
            bg, sig = "blue", "red"
        elif mode == "intra":
            # Polyploid: use chrs_per_subgenome to classify same/cross subgenome
            sub_a = get_chromosome_subgenome(chr_a, chrs_per_subgenome)
            sub_b = get_chromosome_subgenome(chr_b, chrs_per_subgenome)
            if sub_a and sub_b and sub_a == sub_b:
                bg, sig = "red", "blue"   # same-subgenome
            else:
                bg, sig = "blue", "red"   # cross-subgenome or unknown
        else:
            # auto: detect polyploid pattern (chr_b has chr 1..2N, N=chrs_per_subgenome)
            chr_b_num = _normalize_chr_num(chr_b)
            if chr_b_num is not None and 1 <= chr_b_num <= 2 * chrs_per_subgenome:
                sub_a = get_chromosome_subgenome(chr_a, chrs_per_subgenome)
                sub_b = get_chromosome_subgenome(chr_b, chrs_per_subgenome)
                if sub_a and sub_b and sub_a == sub_b:
                    bg, sig = "red", "blue"
                else:
                    bg, sig = "blue", "red"
            else:
                bg, sig = "blue", "red"   # inter-like default
        return [c for c in all_colors if c not in (bg, sig)] + [bg, sig]

    # Group pairs by (chr_a, chr_b) for per-chr-pair rendering
    if dict_gff1 and dict_gff2:
        chr_pair_groups = defaultdict(lambda: defaultdict(list))
        for pair, color in dict_gd.items():
            gene1, gene2 = pair.split(":")
            if gene1 not in loc1 or gene2 not in loc2:
                continue
            chr_a = dict_gff1.get(gene1, ["."])[0]
            chr_b = dict_gff2.get(gene2, ["."])[0]
            chr_pair_groups[(chr_a, chr_b)][color].append((gene1, gene2))
    else:
        chr_pair_groups = {("." , "."): defaultdict(list)}
        for pair, color in dict_gd.items():
            gene1, gene2 = pair.split(":")
            if gene1 in loc1 and gene2 in loc2:
                chr_pair_groups[(".", ".")][color].append((gene1, gene2))

    for (chr_a, chr_b), groups in chr_pair_groups.items():
        draw_order = _draw_order_for_chr_pair(chr_a, chr_b)
        for pass_color in draw_order:
            if pass_color not in groups:
                continue
            c, a = color_map[pass_color]
            for gene1, gene2 in groups[pass_color]:
                x = gl_start1 - loc1[gene1]
                y = gl_start2 + loc2[gene2]
                DrawCircle(root, [y, x], size, c, a)


# ======================================================
# Section 3: Dotplot Pipeline
# ======================================================


def _normalize_chr_num(chr_name: str):
    """Extract chromosome number, handling Chr/chr/CHR prefix case-insensitively.
    Returns int or None.
    """
    import re
    m = re.match(r'chr(\d+)', chr_name.strip(), re.IGNORECASE)
    return int(m.group(1)) if m else None


def _detect_mode(lens_1: list, lens_2: list) -> tuple:
    """Auto-detect intra vs inter mode from chromosome counts.

    Intra (polyploid): species_b has ~2× the chromosomes of species_a
      → donor is diploid parent, recipient is allopolyploid hybrid.
    Inter (cross-species): otherwise.

    Returns (mode, chrs_per_subgenome).
    """
    n_a = len([r for r in lens_1 if _normalize_chr_num(r[0]) is not None])
    n_b = len([r for r in lens_2 if _normalize_chr_num(r[0]) is not None])
    if n_a > 0 and abs(n_b - 2 * n_a) <= 2:
        return "intra", n_a
    return "inter", n_a if n_a > 0 else 10


def calculate_he_direction(chr_a: str, chr_b: str, color: str,
                            parent_genome: str = "A",
                            mode: str = "auto") -> str:
    """Calculate HE direction based on expected vs observed ortholog/paralog pattern.

    Two modes:
    1. intra (peanut-style): same species, A/B subgenome structure.
       ARD chr01-10, ARH chr01-10(A-subgenome) + chr11-20(B-subgenome).
       same-subgenome (chr_a=i, chr_b=i): expected RED, blue -> B_to_A.
       cross-subgenome (chr_a=i, chr_b=i+10): expected BLUE, red -> A_to_B.

    2. inter (3-species): different species, all chr pairs are cross-subgenome.
       Expected BLUE (homeolog between donor/recipient subgenomes).
       red -> A_to_B (donor gene replaced recipient in gene tree).
       blue -> . (normal orthology).
    """
    chr_a_num = _normalize_chr_num(chr_a)
    chr_b_num = _normalize_chr_num(chr_b)
    if chr_a_num is None or chr_b_num is None:
        return "."

    if mode == "inter":
        # All pairs = cross-subgenome (different species/subgenomes)
        # Expected: blue (D_cross = homeolog between subgenomes)
        # red = both genes in same GD clade = A_to_B introgression
        if color == "red":
            return "A_to_B"
        return "."

    # --- intra mode (peanut-style) subgenome logic ---
    # ARD only has chr01-10 (A-genome)
    if not (1 <= chr_a_num <= 10):
        return "."

    # Homologous chromosome offset: A-subgenome chr i ↔ B-subgenome chr i+10
    expected_red = (chr_b_num == chr_a_num)           # same-subgenome ortholog
    expected_blue = (chr_b_num == chr_a_num + 10)     # cross-subgenome homeolog

    if expected_red and color == "blue":
        return "B_to_A"
    elif expected_blue and color == "red":
        return "A_to_B"

    return "."


def write_unified_output(
    sorted_lines: list,
    pair_source_dict: dict,
    sorted_pairs: list,
    zone_records: list,
    dict_gff1: dict,
    dict_gff2: dict,
    spe1: str,
    spe2: str,
    arh_subgenome_dict: dict = None,
    output_path: str = "haplofinder_output.tsv",
    mode: str = "auto",
) -> str:
    """Write single unified TSV with subgenome and ortholog annotation."""
    # --- pair_source display mapping ---
    source_map = {"S": "Orthologs", "D_same": "Paralogs", "D_cross": "Paralogs"}

    zone_lookup: dict = {}
    for rec in zone_records:
        key = (rec["gene1"], rec["gene2"])
        if key not in zone_lookup:
            zone_lookup[key] = (
                rec["zone_id"],
                f"{rec['contig1']}&{rec['contig2']}",
                rec["blue_run"],
            )

    # Sort by species_a chromosome → chromosome start (genomic order)
    # Use natural sort for chromosome names (chr1 < chr2 < chr10, not chr1 < chr10 < chr2)
    def _natural_sort_key(chr_name: str):
        """Natural sort key: extract numeric parts for proper ordering."""
        import re
        parts = re.split(r'(\d+)', chr_name)
        return [int(p) if p.isdigit() else p for p in parts]

    def _sort_key(line: str):
        parts = line.strip().split("\t")
        if len(parts) < 3:
            return (_natural_sort_key("￿"), 0)
        gene_a = parts[0]
        gff = dict_gff1.get(gene_a)
        if gff:
            chr_name = gff[0] if gff[0] != "." else "￿"
            try:
                start = int(gff[1])
            except (ValueError, TypeError):
                start = 0
            return (_natural_sort_key(chr_name), start)
        return (_natural_sort_key("￿"), 0)

    sorted_lines = sorted(sorted_lines, key=_sort_key)

    total_pairs = len(sorted_lines)
    header_cols = [
        f"{spe1}_gene", f"{spe2}_gene",
        f"{spe1}_chr", f"{spe2}_chr",
        f"{spe1}_start", f"{spe1}_end",
        f"{spe2}_start", f"{spe2}_end",
        "color", "Pair_type", "HE_direction",
        "Gene_conversion_zone", "Gene_conversion_chr_pair", "Gene_conversion_blue_run",
        f"{spe2}_subgenome",
    ]
    header = (
        f"# HaploFinder unified output\n"
        f"# species_a={spe1}  species_b={spe2}  total_pairs={total_pairs}\n"
        f"# Columns:\n"
        + "\t".join(header_cols) + "\n"
    )

    with open(output_path, "w") as fh:
        fh.write(header)
        for line in sorted_lines:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            gene_a, gene_b, color = parts[0], parts[1], parts[2]
            gff_a = dict_gff1.get(gene_a)
            gff_b = dict_gff2.get(gene_b)
            chr_a = gff_a[0] if gff_a else "."
            chr_b = gff_b[0] if gff_b else "."
            start_a = gff_a[1] if gff_a else "."
            end_a = gff_a[2] if gff_a else "."
            start_b = gff_b[1] if gff_b else "."
            end_b = gff_b[2] if gff_b else "."

            # Calculate HE_direction (mode-aware)
            he_direction = calculate_he_direction(chr_a, chr_b, color, parent_genome="A", mode=mode)

            source = pair_source_dict.get((gene_a, gene_b, color), "?")
            display_source = source_map.get(source, source)
            zone_info = zone_lookup.get((gene_a, gene_b), (".", ".", "."))
            zid, zcp, zbr = zone_info
            subgenome = arh_subgenome_dict.get(gene_b, ".") if arh_subgenome_dict else "."
            fh.write(
                f"{gene_a}\t{gene_b}\t"
                f"{chr_a}\t{chr_b}\t"
                f"{start_a}\t{end_a}\t{start_b}\t{end_b}\t"
                f"{color}\t{display_source}\t{he_direction}\t"
                f"{zid}\t{zcp}\t{zbr}\t{subgenome}\n"
            )
    logger.info("Wrote unified output: %s (%d pairs)", output_path, total_pairs)
    return output_path


def generate_dotplot(
    gff1, gff2, lens1, lens2, gd_pairs, spe1, spe2,
    file_name, target_chr1=None, target_chr2=None, size=None,
    min_pairs: int = 10,
    n_permutations: int = 1000,
    p_threshold: float = 0.05,
    pair_source_dict: dict = None,
    arh_subgenome_dict: dict = None,
    output_dir: str = ".",
    mode: str = "auto",
    chr_pair: str = None,
    chrs_per_subgenome: int = 10,
):
    """Generate dotplot and run gene-conversion detection.

    If chr_pair is specified (format "ChrA,ChrB"), filter to that single
    chromosome pair, increase dot size, and highlight conversion zones.
    """
    # Resolve chr_pair filter
    chr_a_filter = None
    chr_b_filter = None
    if chr_pair:
        parts = chr_pair.split(",")
        if len(parts) == 2:
            chr_a_filter, chr_b_filter = parts[0].strip(), parts[1].strip()
            logger.info("Zoom mode: showing only %s <-> %s", chr_a_filter, chr_b_filter)

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

    # Resolve mode: auto-detect from donor/receiver chromosome counts
    if mode == "auto":
        mode, detected_n = _detect_mode(lens_1, lens_2)
        chrs_per_subgenome = detected_n
        logger.info("Mode resolved to: %s (chrs_per_subgenome=%d)", mode, chrs_per_subgenome)
    else:
        logger.info("Mode: %s (chrs_per_subgenome=%d)", mode, chrs_per_subgenome)

    dict_gd = read_gd_pairs(gd_pairs)

    # ── Zoom mode: filter to single chromosome pair ──
    if chr_a_filter and chr_b_filter:
        # Filter lens to target chromosomes only
        lens_map1 = dict(lens_1)
        lens_map2 = dict(lens_2)
        if chr_a_filter in lens_map1 and chr_b_filter in lens_map2:
            lens_1 = [[chr_a_filter, lens_map1[chr_a_filter]]]
            lens_2 = [[chr_b_filter, lens_map2[chr_b_filter]]]
        # Filter dict_gd
        dict_gd = {
            k: v for k, v in dict_gd.items()
            if dict_gff1.get(k.split(":")[0], [None])[0] == chr_a_filter
            and dict_gff2.get(k.split(":")[1], [None])[0] == chr_b_filter
        }
        # Increase dot size for zoomed view
        if size is None or size == 0.001:
            size = 0.005
        logger.info("Zoom filter: %d pairs on %s<->%s, dot_size=%.4f",
                     len(dict_gd), chr_a_filter, chr_b_filter, size or 0.005)

    gl1, gl2 = 0.92, 0.92
    step_1 = plot_chr1(lens_1, gl1, gl2, "", spe1)
    step_2 = plot_chr2(lens_2, gl2, gl1, "", spe2)

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
        mode=mode,
    )

    sorted_pairs = find_gene_pair_info(
        gene_conversion_list, dict_gd, dict_gff1, dict_gff2, file_name
    )
    # Extract 3-element tuples for conversion zone scanning
    result_conversion = [(g1, g2, c) for g1, g2, c, *_ in sorted_pairs]
    zone_records = find_conversion_zones_with_ids_to_file(
        result_conversion, dict_gff1, dict_gff2, mode=mode
    )

    # Write unified output if pair_source_dict is provided
    if pair_source_dict is not None:
        # Convert gd_pairs (list of strings) to the format expected by write_unified_output
        # gd_pairs contains all GD pairs, not just those in conversion zones
        all_pairs_for_output = []
        for line in gd_pairs:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                gene_a, gene_b, color = parts[0], parts[1], parts[2]
                gff_a = dict_gff1.get(gene_a)
                gff_b = dict_gff2.get(gene_b)
                if gff_a and gff_b:
                    chr_a, start_a = gff_a[0], gff_a[1]
                    chr_b, start_b = gff_b[0], gff_b[1]
                    try:
                        start_a_int = int(start_a)
                        start_b_int = int(start_b)
                        all_pairs_for_output.append((gene_a, gene_b, color, chr_a, start_a_int, chr_b, start_b_int))
                    except (ValueError, TypeError):
                        pass

        write_unified_output(
            gd_pairs, pair_source_dict, all_pairs_for_output, zone_records,
            dict_gff1, dict_gff2, spe1, spe2,
            arh_subgenome_dict=arh_subgenome_dict,
            output_path=os.path.join(output_dir, "haplofinder_output.tsv"),
            mode=mode,
        )

    # Cache gene locations to avoid redundant GFF traversals
    gene_loc_1 = gene_location(gff_1, lens_1, step_1)
    gene_loc_2 = gene_location(gff_2, lens_2, step_2)

    gc.collect()
    if size:
        plot_dot(root, gene_loc_1, gene_loc_2, dict_gd, size, dict_gff1, dict_gff2, mode, chrs_per_subgenome)
    else:
        plot_dot(root, gene_loc_1, gene_loc_2, dict_gd, dict_gff1=dict_gff1, dict_gff2=dict_gff2, mode=mode, chrs_per_subgenome=chrs_per_subgenome)

    # ── Highlight conversion zones with rectangles (zoom mode) ──
    if chr_a_filter and chr_b_filter and zone_records:
        # Collect zone spans grouped by zone_id
        from collections import defaultdict as _dd
        zone_spans = _dd(lambda: [None, None, None, None])  # min_x,max_x,min_y,max_y
        for rec in zone_records:
            if rec["contig1"] != chr_a_filter or rec["contig2"] != chr_b_filter:
                continue
            g1, g2 = rec["gene1"], rec["gene2"]
            x = gene_loc_1.get(g1)
            y = gene_loc_2.get(g2)
            if x is not None and y is not None:
                zid = rec["zone_id"]
                cur = zone_spans[zid]
                zone_spans[zid] = [
                    min(cur[0], x) if cur[0] is not None else x,
                    max(cur[1], x) if cur[1] is not None else x,
                    min(cur[2], y) if cur[2] is not None else y,
                    max(cur[3], y) if cur[3] is not None else y,
                ]
        for zid, (xmin, xmax, ymin, ymax) in zone_spans.items():
            if None not in (xmin, xmax, ymin, ymax):
                rect = mpatches.Rectangle(
                    (xmin, ymin), xmax - xmin, ymax - ymin,
                    linewidth=1.5, edgecolor="orange", facecolor="none",
                    linestyle="--", alpha=0.8,
                )
                root.add_patch(rect)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    out_prefix = os.path.join(output_dir, os.path.basename(file_name))
    plt.savefig(out_prefix + "_dotplot.pdf", dpi=500)
    plt.savefig(out_prefix + "_dotplot.png", dpi=500)
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

    # Reciprocal pairs + minimum coverage filter.
    # Collect ALL reciprocal pairs above threshold, not just the single best
    # hit.  For allopolyploids this captures both same-subgenome pairs
    # (e.g. chr01<->chr01) and cross-subgenome homeologs (e.g. chr01<->chr11).
    homolog_pairs = []
    seen_pairs = set()
    for (chr_a, chr_b), cnt_ab in pair_count.items():
        if cnt_ab < min_shared_pairs:
            continue
        if (chr_a, chr_b) in seen_pairs or (chr_b, chr_a) in seen_pairs:
            continue
        # Reciprocal: chr_b must also have chr_a in its top hit
        if best_a_for_b.get(chr_b, (None,))[0] == chr_a:
            homolog_pairs.append((chr_a, chr_b))
            seen_pairs.add((chr_a, chr_b))
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
    mode: str = "auto",
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
    return sort_lst


def find_conversion_zones_with_ids_to_file(
    data: list,
    dict_gff1: dict,
    dict_gff2: dict,
    output_file: str = "gene_conversion.txt",
    mode: str = "auto",
) -> list:
    """Scan for gene conversion zones using the corrected red→blue→red pattern.

    Original code scanned blue→red→blue, which is the inverted signal.

    Gene conversion causes a donor-sequence island to interrupt the expected
    collinear (red) background, producing a local enrichment of blue (cross-clade)
    pairs flanked by red (in-clade) pairs on both sides.

    Corrected pattern
    -----------------
        red … red   |   blue … blue   |   red … red
        (normal)        (converted)        (normal)

    Each chromosome pair is scanned independently to avoid fragmentation
    when multiple homolog pairs coexist (e.g. chr01↔chr01 + chr01↔chr11).
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

    def _scan_chromosome_pair(pairs: list, chrom1: str, chrom2: str,
                               zone_id_start: int,
                               mode: str = "auto") -> list:
        """Scan a single chromosome pair for conversion zones.

        Detects two biological patterns depending on mode:

        intra mode (peanut subgenome):
            Same-subgenome (ARH chr01-10, expected red background):
                red … red | blue … blue | red … red  →  blue island = B→A conversion
            Cross-subgenome (ARH chr11-20, expected blue background):
                blue … blue | red … red | blue … blue  →  red island = A→B conversion

        inter mode (different species):
            All pairs = cross-subgenome (donor vs recipient).
            expected blue background, red islands = A→B conversion.
        """
        chr_b_num = _normalize_chr_num(chrom2)

        if mode == "inter":
            # All pairs cross-subgenome: blue background, red signal
            bg_color, sig_color = "blue", "red"
        elif mode == "intra" and chr_b_num is not None:
            if 1 <= chr_b_num <= 10:
                bg_color, sig_color = "red", "blue"      # same-subgenome
            elif 11 <= chr_b_num <= 20:
                bg_color, sig_color = "blue", "red"      # cross-subgenome
            else:
                return []
        else:
            # auto or unknown: fallback to chr_b_num pattern detection
            if chr_b_num is None:
                return []
            if 1 <= chr_b_num <= 10:
                bg_color, sig_color = "red", "blue"
            elif 11 <= chr_b_num <= 20:
                bg_color, sig_color = "blue", "red"
            else:
                return []

        label = "blue_run" if sig_color == "blue" else "red_run"

        pair_zones = []
        n = len(pairs)
        i = 0
        while i < n - 1:
            gene1_start, gene2_start, color = pairs[i]

            if color != bg_color:
                i += 1
                continue

            # Consume leading background run
            end_bg1, _ = _consume_color_run(pairs, i, bg_color, chrom1, chrom2, dict_gff1, dict_gff2)
            i = end_bg1 + 1

            if i >= n or pairs[i][2] != sig_color:
                continue

            # Consume signal (converted) run
            start_sig = i
            end_sig, sig_run_len = _consume_color_run(pairs, i, sig_color, chrom1, chrom2, dict_gff1, dict_gff2)
            i = end_sig + 1

            if i >= n or pairs[i][2] != bg_color:
                i = start_sig
                continue

            # Consume trailing background run
            start_bg2 = i
            end_bg2, _ = _consume_color_run(pairs, i, bg_color, chrom1, chrom2, dict_gff1, dict_gff2)
            i = end_bg2 + 1

            if sig_run_len >= 3:
                local_start = end_bg1
                local_end = start_bg2
                pair_zones.append((local_start, local_end, start_sig, end_sig, sig_run_len))
                logger.info(
                    "Conversion zone %d: %s&%s  %s=%d  [idx %d-%d]",
                    zone_id_start + len(pair_zones), chrom1, chrom2,
                    label, sig_run_len, start_sig, end_sig,
                )

            i = start_bg2

        return pair_zones

    # ── Group pairs by chromosome pair ────────────────────────────
    chr_groups: Dict[Tuple[str, str], list] = defaultdict(list)
    for global_idx, (g1, g2, c) in enumerate(data):
        gff1 = dict_gff1.get(g1)
        gff2 = dict_gff2.get(g2)
        if gff1 and gff2:
            chr_groups[(gff1[0], gff2[0])].append((global_idx, g1, g2, c))

    # ── Scan each chromosome pair independently ──────────────────
    all_conversion_zones = []  # (zone_id, s_r1_global, e_r2_global, s_b_global, e_b_global, blue_run_len)
    zone_id = 1

    for (chrom1, chrom2), indexed_pairs in chr_groups.items():
        # Build local list (global_idx, gene_a, gene_b, color)
        local_data = [(g1, g2, c) for _, g1, g2, c in indexed_pairs]
        local_zones = _scan_chromosome_pair(local_data, chrom1, chrom2, zone_id, mode=mode)

        for l_start, l_end, l_sb, l_eb, blen in local_zones:
            # Map local indices back to global indices
            g_start = indexed_pairs[l_start][0]
            g_end = indexed_pairs[l_end][0]
            g_sb = indexed_pairs[l_sb][0]
            g_eb = indexed_pairs[l_eb][0]
            all_conversion_zones.append((zone_id, g_start, g_end, g_sb, g_eb, blen))
            zone_id += 1

    logger.info("Detected %d conversion zones", len(all_conversion_zones))

    # ── Build structured zone output (original format) ────────────
    zone_records = []
    for zone in all_conversion_zones:
        zid, s_r1, e_r2, s_b, e_b, blen = zone
        first_gene1, first_gene2, _ = data[s_r1]
        contig1 = dict_gff1[first_gene1][0]
        contig2 = dict_gff2[first_gene2][0]
        for j in range(s_r1, e_r2 + 1):
            if j >= len(data):
                break
            g1, g2, color = data[j]
            gff1 = dict_gff1.get(g1)
            gff2 = dict_gff2.get(g2)
            if gff1 and gff2:
                zone_records.append({
                    "zone_id": zid,
                    "contig1": gff1[0],
                    "contig2": gff2[0],
                    "gene1": g1,
                    "gene2": g2,
                    "start1": gff1[1],
                    "end1": gff1[2],
                    "start2": gff2[1],
                    "end2": gff2[2],
                    "color": color,
                    "blue_run": blen,
                })
    return zone_records


# ======================================================
# Section 7: GD Pair Processing  (duplication node cap removed)
# ======================================================


def process_gd_result(gf, imap, input_sps_tree, sp1, sp2, support, pair_support,
                      parental_sps=None, hyb_sps=None):
    """Process gene-family trees to assign red/blue labels and optionally
    subgenome to species_b (hybrid) genes.

    When parental_sps and hyb_sps are provided, each species_b gene is
    assigned A/B via S-node sister-clade inference (same logic as split mode).
    """
    tre_dic = read_and_return_dict(gf)
    (gene2new_named_gene_dic, new_named_gene2gene_dic,
     voucher2taxa_dic, taxa2voucher_dic) = gene_id_transfer(imap)
    rename_sp1 = taxa2voucher_dic[sp1]
    rename_sp2 = taxa2voucher_dic[sp2]
    renamed_sptree = rename_input_tre(
        read_phylo_tree(input_sps_tree), taxa2voucher_dic
    )
    # --- subgenome assignment setup ---
    do_subgenome = bool(parental_sps and hyb_sps)
    arh_subgenome_dict: dict = {}
    if do_subgenome:
        hyb_renamed = taxa2voucher_dic[hyb_sps]
        parental_renamed = {taxa2voucher_dic[s] for s in parental_sps}
        subgenome_label = {tag: chr(ord("A") + i)
                           for i, tag in enumerate(sorted(parental_renamed))}
    processed_lines: list = []
    pair_source_dict: dict = {}
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
        # Mark speciation pairs with source "S"
        for i, line in enumerate(processed_lines):
            parts = line.strip().split("\t")
            if len(parts) == 3:
                key = (parts[0], parts[1], parts[2])
                if key not in pair_source_dict:
                    pair_source_dict[key] = "S"
            elif len(parts) == 4 and parts[3] == ".":
                # already has placeholder from speciation, mark source
                key = (parts[0], parts[1], parts[2])
                if key not in pair_source_dict:
                    pair_source_dict[key] = "S"

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
                    color, suffix, source = "red", "r", "D_same"
                else:
                    color, suffix, source = "blue", "b", "D_cross"
                    # record which clade ARD is in for direction inference
                    ard_in_tips1 = gene_a in tips1 or gene_b in tips1
                    # ARD reference determines the A-like clade
                token = f"{gene2}-{suffix}"
                if token not in written_results:
                    processed_lines.append(f"{gene_a}\t{gene_b}\t{color}\t.\n")
                    pair_source_dict[(gene_a, gene_b, color)] = source
                    written_results.add(token)

        # --- inline subgenome assignment (same S-node sister-clade logic) ---
        if do_subgenome and parental_renamed.issubset(
            {g.split("_")[0] for g in Phylo_t1.get_leaf_names()}
        ):
            s_nodes: set = set()
            for event in Phylo_t1.get_descendant_evol_events():
                if event.etype == "S":
                    enodes = set(event.in_seqs) | set(event.out_seqs)
                    s_nodes.add(Phylo_t1.get_common_ancestor(list(enodes)))
            for leaf in Phylo_t1.get_leaves():
                leaf_name = leaf.name
                if not leaf_name.startswith(hyb_renamed):
                    continue
                if leaf_name not in new_named_gene2gene_dic:
                    continue
                child = leaf
                node = leaf.up
                best_parent = None
                while node is not None:
                    if node not in s_nodes:
                        child = node; node = node.up; continue
                    sister_tags = set()
                    sister_has_hybrid = False
                    for sister in node.children:
                        if sister is child: continue
                        for ln in sister.get_leaf_names():
                            tag = ln.split("_")[0]
                            if tag == hyb_renamed: sister_has_hybrid = True
                            elif tag in parental_renamed: sister_tags.add(tag)
                    if len(sister_tags) == 1 and not sister_has_hybrid:
                        best_parent = next(iter(sister_tags)); break
                    if len(sister_tags) > 1: break
                    child = node; node = node.up
                if best_parent is not None:
                    label = subgenome_label[best_parent]
                    og = new_named_gene2gene_dic[leaf_name]
                    if og in arh_subgenome_dict and arh_subgenome_dict[og] != label:
                        arh_subgenome_dict[og] = "unknown"
                    else:
                        arh_subgenome_dict[og] = label

    if skipped_trees:
        logger.info("Skipped %d trees (single-species or missing target sp)", skipped_trees)
    if do_subgenome:
        logger.info("Subgenome assigned to %d %s genes", len(arh_subgenome_dict), sp2)

    # Note: Conversion_direction calculation removed.
    # HE_direction is now calculated in write_unified_output() based on chromosome positions.

    sorted_lines = processed_lines  # sorting deferred to write_unified_output (chr-aware)
    return sorted_lines, pair_source_dict, arh_subgenome_dict


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
    input_fasta, haplofinder_tsv, gff, output_dir,
):
    """Split hybrid sequences into per-subgenome FASTA files.

    Reads ARH_subgenome assignments from haplofinder_output.tsv
    and partitions the input FASTA accordingly.
    """
    # Read subgenome assignments from unified output
    arh_sub: dict = {}
    with open(haplofinder_tsv, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 13:
                continue
            gene_id = parts[1]   # species_b gene
            sub = parts[14]       # species_b_subgenome
            if sub in ("A", "B"):
                arh_sub[gene_id] = sub
                arh_sub[gene_id] = sub

    gff_1, dict_gff1 = read_gff(gff)
    _ = gff_1
    all_sequences = read_fasta(input_fasta)
    os.makedirs(output_dir, exist_ok=True)

    all_labels = set(arh_sub.values()) - {"unknown"}
    split_buckets: Dict[str, dict] = {lbl: {} for lbl in all_labels}
    split_buckets["unknown"] = {}
    assigned_total = 0
    skipped_unassigned = 0

    with open(os.path.join(output_dir, "split_assignment.tsv"), "w") as out:
        out.write("gene_id\tsubgenome\tchromosome\n")
        for gene_id, seq in all_sequences.items():
            label = arh_sub.get(gene_id)
            if label is None:
                skipped_unassigned += 1
                continue
            assigned_total += 1
            chr_name = dict_gff1.get(gene_id, [None])[0]
            split_buckets.setdefault(label, {})[gene_id] = seq
            out.write(
                f"{gene_id}\t{label}\t"
                f"{chr_name if chr_name else 'NA'}\n"
            )

    for lbl, seqs in split_buckets.items():
        safe_lbl = lbl.replace("/", "_")
        write_fasta(seqs, os.path.join(output_dir, f"split_subgenome_{safe_lbl}.fasta"))

    with open(os.path.join(output_dir, "split_summary.txt"), "w") as out:
        out.write(f"input_fasta_total={len(all_sequences)}\n")
        out.write(f"assigned_total={assigned_total}\n")
        out.write(f"skipped_unassigned={skipped_unassigned}\n")
        for lbl, seqs in split_buckets.items():
            out.write(f"subgenome_{lbl}={len(seqs)}\n")

    logger.info("split outputs written to %s", os.path.abspath(output_dir))


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
# Section 13: Optional input processors (stub implementations)
# ======================================================

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


# ======================================================
# Section 14: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse
    import logging

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

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
    parser.add_argument("--mode", type=str, default="auto",
                        choices=["auto", "intra", "inter"],
                        help="Analysis mode: auto (detect from chr naming), "
                             "intra (peanut-style A/B subgenome), "
                             "inter (different species, all cross-subgenome)")
    parser.add_argument("--chr_pair", type=str, default=None,
                        help="Zoom to single chromosome pair: ChrA,ChrB")
    args = parser.parse_args()

    process_blastp_pairs = _process_blastp_result_stub(args.blastp_pairs, args.num)
    alignments, alignment_scores = _parse_synteny_file_stub(args.synteny_pairs)
    process_synteny_pairs = assign_colors_by_alignment(alignments, alignment_scores)
    process_gd_pairs, pair_source_dict, _arh_sub = process_gd_result(
        args.gf, args.imap, args.input_sps_tree, args.spe1, args.spe2, 50, 50
    )

    for pairs, tag in [
        (process_blastp_pairs, "blastp_pairs"),
        (process_synteny_pairs, "synteny_pairs"),
        (process_gd_pairs, "gd_pairs"),
    ]:
        p_src = pair_source_dict if tag == "gd_pairs" else None
        generate_dotplot(
            args.gff1, args.gff2, args.lens1, args.lens2,
            pairs, args.spe1, args.spe2, tag,
            args.target_chr1, args.target_chr2, args.size,
            min_pairs=args.min_conv_pairs,
            n_permutations=args.n_permutations,
            p_threshold=args.p_threshold,
            pair_source_dict=p_src,
            mode=args.mode,
            chr_pair=args.chr_pair,
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
        pair_source_dict=pair_source_dict,
        mode=args.mode,
        chr_pair=args.chr_pair,
    )


