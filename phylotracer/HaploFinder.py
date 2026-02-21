#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Haplotype detection and gene conversion analysis for the PhyloTracer pipeline.

This module supports dotplot visualization, duplication-based gene labeling,
conversion zone detection, and hybrid subgenome assignment.
"""

import gc
import importlib
import os
import re
import sys
import time
from collections import Counter, defaultdict

import matplotlib
import matplotlib.patches as mpatches
from ete3 import PhyloTree
from tqdm import tqdm

from phylotracer import (
    annotate_gene_tree,
    find_dup_node,
    gene_id_transfer,
    get_species_set,
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

# ======================================================
# Section 1: GFF and Lens Parsing
# ======================================================


def read_gff(fn):
    """
    Read a simplified GFF-like file into a list and dictionary.

    Parameters
    ----------
    fn : str
        Path to a GFF-like input file.

    Returns
    -------
    tuple
        (data_list, gene_dict) where gene_dict maps gene IDs to attributes.

    Assumptions
    -----------
    Input rows are tab-delimited with fields including gene identifiers.
    """
    with open(fn) as f:
        data, data_dict = [], {}
        for line in f.readlines():
            a = line.strip().split("\t")
            data_dict[a[1]] = [a[0], a[2], a[3], a[4], a[5]]
            data.append(a)
    return data, data_dict


def read_lens(fn, chrs=None):
    """
    Read chromosome length data with optional filtering.

    Parameters
    ----------
    fn : str
        Path to a lens file.
    chrs : str, optional
        Path to a list of chromosome names to retain.

    Returns
    -------
    list
        List of [chromosome, length] entries.

    Assumptions
    -----------
    Lens file is tab-delimited with chromosome and length columns.
    """
    data = []

    chrs_lst = None
    if chrs:
        with open(chrs, "r") as f:
            chrs_lst = {i.strip() for i in f.readlines()}

    with open(fn, "r") as fp:
        for row in fp.readlines():
            r1, r2 = row.strip().split("\t")
            if chrs_lst is None or r1 in chrs_lst:
                data.append([r1, r2])

    return data


# ======================================================
# Section 2: Dotplot Rendering Helpers
# ======================================================


def plot_chr1(lens, gl, gl2, mark, name):
    """
    Plot chromosome segments for the first species axis.

    Parameters
    ----------
    lens : list
        Chromosome length entries.
    gl : float
        Total plot length for the axis.
    gl2 : float
        Secondary length for drawing.
    mark : object
        Prefix for chromosome labels.
    name : str
        Species name for axis labeling.

    Returns
    -------
    float
        Step size for coordinate mapping.

    Assumptions
    -----------
    Matplotlib axes are already initialized.
    """
    total_lens = sum([float(k[1]) for k in lens])
    step = gl / float(total_lens)
    gl_start, n, start_x = 0.95, 0, 0.05
    mark_y = 0.04
    align = dict(family="Times New Roman", style="normal", horizontalalignment="center", verticalalignment="center")
    for k in lens:
        n += float(k[1])
        mark_new = str(mark) + str(k[0])
        x = gl_start - float(n) * step
        mark_x = x + 0.5 * float(k[1]) * step
        plt.plot([start_x, start_x + gl2], [x, x], linestyle="-", color="black", linewidth=0.5)
        plt.text(mark_y, mark_x, mark_new, color="black", fontsize=12, rotation=90, weight="semibold", **align)
    plt.plot([start_x, start_x + gl2], [gl_start, gl_start], linestyle="-", color="black", linewidth=1)
    plt.text(mark_y - 0.02, 0.5 * (2 * gl_start - gl), name, color="black", fontsize=18, rotation=90, weight="semibold", **align)
    return step


def plot_chr2(lens, gl, gl2, mark, name):
    """
    Plot chromosome segments for the second species axis.

    Parameters
    ----------
    lens : list
        Chromosome length entries.
    gl : float
        Total plot length for the axis.
    gl2 : float
        Secondary length for drawing.
    mark : object
        Prefix for chromosome labels.
    name : str
        Species name for axis labeling.

    Returns
    -------
    float
        Step size for coordinate mapping.

    Assumptions
    -----------
    Matplotlib axes are already initialized.
    """
    total_lens = sum([float(k[1]) for k in lens])
    step = gl / float(total_lens)
    gl_start, n, start_x = 0.05, 0, 0.95
    mark_y = 0.96
    align = dict(family="Times New Roman", style="normal", horizontalalignment="center", verticalalignment="center")
    for k in lens:
        n += float(k[1])
        mark_new = str(mark) + str(k[0])
        x = gl_start + float(n) * step
        mark_x = x - 0.5 * float(k[1]) * step
        plt.plot([x, x], [start_x, start_x - gl2], linestyle="-", color="black", linewidth=0.5)
        plt.text(mark_x, mark_y, mark_new, color="black", fontsize=12, rotation=0, weight="semibold", **align)
    plt.plot([gl_start, gl_start], [start_x, start_x - gl2], linestyle="-", color="black", linewidth=1)
    plt.text(0.5 * (2 * gl_start + gl), mark_y + 0.02, name, color="black", fontsize=18, rotation=0, weight="semibold", **align)
    return step


def gene_location(gff, lens, step):
    """
    Map genes to scaled genomic coordinates.

    Parameters
    ----------
    gff : list
        GFF-like entries.
    lens : list
        Chromosome length entries.
    step : float
        Scaling factor for coordinates.

    Returns
    -------
    dict
        Mapping from gene ID to scaled location.

    Assumptions
    -----------
    Gene entries store chromosome and positional information.
    """
    loc_gene, dict_chr, n = {}, {}, 0
    for i in lens:
        dict_chr[i[0]] = n
        n += float(i[1])
    for k in gff:
        if k[0] not in dict_chr.keys():
            continue
        loc = (float(dict_chr[k[0]]) + float(k[5])) * step
        loc_gene[k[1]] = loc
    return loc_gene


def read_gd_pairs(gd_lst):
    """
    Parse GD pair labels into a symmetric dictionary.

    Parameters
    ----------
    gd_lst : list
        List of GD pair strings.

    Returns
    -------
    dict
        Mapping from gene pair keys to color labels.

    Assumptions
    -----------
    Input lines are tab-delimited with gene1, gene2, and label.
    """
    dict_gd, dict_gd1 = {}, {}
    for line in gd_lst:
        a = line.strip().split("\t")
        dict_gd1[str(a[0]) + ":" + str(a[1])] = str(a[2])
    for pair, wgd_notch in dict_gd1.items():
        fd = pair.strip().split(":")
        gene1 = fd[0]
        gene2 = fd[1]
        index = gene1 + ":" + gene2
        reverse_index = gene2 + ":" + gene1
        dict_gd[index] = dict_gd1[index]
        if reverse_index not in dict_gd1.keys():
            dict_gd[reverse_index] = dict_gd1[index]
    return dict_gd


def plot_dot(root, loc1, loc2, dict_gd, size=0.001):
    """
    Plot GD pairs as colored dots in a dotplot.

    Parameters
    ----------
    root : object
        Matplotlib axes for plotting.
    loc1 : dict
        Gene locations for species 1.
    loc2 : dict
        Gene locations for species 2.
    dict_gd : dict
        Mapping from gene pair keys to color labels.
    size : float, optional
        Dot radius.

    Returns
    -------
    None

    Assumptions
    -----------
    Locations are scaled between 0 and 1.
    """
    gl_start1, gl_start2 = 0.95, 0.05

    for pair, wgd_notch in dict_gd.items():
        fd = pair.strip().split(":")
        gene1 = fd[0]
        gene2 = fd[1]
        index = gene1 + ":" + gene2
        if gene1 in loc1.keys() and gene2 in loc2.keys():
            x, y = loc1[gene1], loc2[gene2]
            x, y = gl_start1 - x, gl_start2 + y
            if dict_gd[index] == "NONE":
                DrawCircle(root, [y, x], size, "gray", 0.6)

    for pair, wgd_notch in dict_gd.items():
        fd = pair.strip().split(":")
        gene1 = fd[0]
        gene2 = fd[1]
        index = gene1 + ":" + gene2
        if gene1 in loc1.keys() and gene2 in loc2.keys():
            x, y = loc1[gene1], loc2[gene2]
            x, y = gl_start1 - x, gl_start2 + y
            if dict_gd[index] == "red":
                DrawCircle(root, [y, x], size, "red", 0.6)
            elif dict_gd[index] == "green":
                DrawCircle(root, [y, x], size, "green", 0.6)
            elif dict_gd[index] == "purple":
                DrawCircle(root, [y, x], size, "purple", 0.6)
            elif dict_gd[index] == "blue":
                DrawCircle(root, [y, x], size, "blue", 0.6)


def DrawCircle(ax, loc, radius, color, alpha):
    """
    Draw a filled circle on a matplotlib axes.

    Parameters
    ----------
    ax : object
        Matplotlib axes.
    loc : list
        [x, y] coordinate of the circle center.
    radius : float
        Circle radius.
    color : str
        Circle color.
    alpha : float
        Transparency value.

    Returns
    -------
    None

    Assumptions
    -----------
    Axes supports patch addition.
    """
    circle = mpatches.Circle(loc, radius, edgecolor="none", facecolor=color, alpha=alpha)
    ax.add_patch(circle)


# ======================================================
# Section 3: Dotplot Pipeline
# ======================================================


def generate_dotplot(
    gff1,
    gff2,
    lens1,
    lens2,
    gd_pairs,
    spe1,
    spe2,
    file_name,
    target_chr1=None,
    target_chr2=None,
    size=None,
):
    """
    Generate a dotplot for gene pairs with optional chromosome filtering.

    Parameters
    ----------
    gff1 : str
        GFF-like file for species 1.
    gff2 : str
        GFF-like file for species 2.
    lens1 : str
        Lens file for species 1.
    lens2 : str
        Lens file for species 2.
    gd_pairs : list
        GD pair list or label list.
    spe1 : str
        Species 1 label.
    spe2 : str
        Species 2 label.
    file_name : str
        Output file prefix.
    target_chr1 : str, optional
        Chromosome filter for species 1.
    target_chr2 : str, optional
        Chromosome filter for species 2.
    size : float, optional
        Dot size override.

    Returns
    -------
    None

    Assumptions
    -----------
    Input files are formatted according to the expected dotplot schema.
    """
    plt.figure(figsize=(10, 10))
    root = plt.axes([0, 0, 1, 1])
    align = dict(family="Arial", style="normal", horizontalalignment="center", verticalalignment="center")
    t1 = time.time()
    print(f"Dotplot of {file_name} are ready to begin")
    gff_1, dict_gff1 = read_gff(gff1)
    gff_2, dict_gff2 = read_gff(gff2)
    t2 = time.time()
    print("Reading gff took " + str(t2 - t1) + " second")
    if target_chr1 and target_chr2:
        lens_1 = read_lens(lens1, target_chr1)
        lens_2 = read_lens(lens2, target_chr2)
    else:
        lens_1 = read_lens(lens1)
        lens_2 = read_lens(lens2)
    t3 = time.time()
    print("Reading lens took " + str(t3 - t2) + " second")
    gl1, gl2 = 0.92, 0.92
    step_1 = plot_chr1(lens_1, gl1, gl2, "", spe1)
    step_2 = plot_chr2(lens_2, gl2, gl1, "", spe2)

    dict_gd = read_gd_pairs(gd_pairs)
    t4 = time.time()
    print("Reading lebel_pairs took " + str(t4 - t3) + " second")
    gene_loc_1 = gene_location(gff_1, lens_1, step_1)
    gene_loc_2 = gene_location(gff_2, lens_2, step_2)

    gene_conversion_list = find_gene_conversion(
        dict_gd,
        dict_gff1,
        dict_gff2,
        lens_1,
        lens_2,
    )

    result_conversion = find_gene_pair_info(gene_conversion_list, dict_gd, dict_gff1, dict_gff2, file_name)
    find_conversion_zones_with_ids_to_file(result_conversion, dict_gff1, dict_gff2)
    t5 = time.time()
    print("Dealing lebel_pairs took " + str(t5 - t4) + " second")
    gc.collect()

    if size:
        plot_dot(root, gene_loc_1, gene_loc_2, dict_gd, size)
    else:
        plot_dot(root, gene_loc_1, gene_loc_2, dict_gd)
    t6 = time.time()
    print("Ploting dot took " + str(t6 - t5) + " second")
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    plt.savefig(file_name + "_dotplot.pdf", dpi=500)
    plt.savefig(file_name + "_dotplot.png", dpi=500)
    t7 = time.time()
    print(f"{file_name} Dotplot totaly took " + str(t7 - t1) + " second")


# ======================================================
# Section 4: Alignment Labeling
# ======================================================


def assign_colors_by_alignment(alignments, alignment_scores):
    """
    Assign colors to gene pairs based on alignment ranking.

    Parameters
    ----------
    alignments : dict
        Mapping from alignment identifiers to gene-pair lists.
    alignment_scores : dict
        Mapping from alignment identifiers to scores.

    Returns
    -------
    list
        Sorted list of color-labeled gene pairs.

    Assumptions
    -----------
    Higher scores indicate stronger alignments.
    """
    gene_to_alignments = defaultdict(list)

    for alignment, gene_pairs in alignments.items():
        for gene_a, gene_b in gene_pairs:
            gene_to_alignments[gene_a].append((alignment, gene_b))

    result = set()
    gene_pair_set = set()

    for gene_a, alignments_list in gene_to_alignments.items():
        alignments_list.sort(key=lambda x: alignment_scores[x[0]], reverse=True)

        for i, (alignment, gene_b) in enumerate(alignments_list):
            if i == 0:
                color = "red"
            else:
                color = "blue"
            gene_pair = f"{gene_a}\t{gene_b}"
            if gene_pair not in gene_pair_set:
                result.add(f"{gene_pair}\t{color}")
                gene_pair_set.add(gene_pair)
            else:
                continue

    sorted_list = sorted(result)
    return sorted_list


def judge_support(support, support_value):
    """
    Evaluate whether a support value meets a threshold.

    Parameters
    ----------
    support : float
        Support value to evaluate.
    support_value : float
        Threshold value.

    Returns
    -------
    bool
        True if the support passes the threshold.

    Assumptions
    -----------
    Support values may be expressed as fractions or percentages.
    """
    if support <= 1 and 0.5 <= support_value <= 1:
        if support >= support_value:
            return True
        else:
            return False

    elif support <= 1 and 50 <= support_value <= 100:
        support_value = support_value / 100
        if support >= support_value:
            return True
        else:
            return False
    elif support > 1 and 0.5 <= support_value <= 1:
        support_value = support_value * 100
        if support >= support_value:
            return True
        else:
            return False
    elif support > 1 and 50 <= support_value <= 100:
        if support >= support_value:
            return True
        else:
            return False


# ======================================================
# Section 5: Duplication and Orthology Detection
# ======================================================


def generate_combinations(list1, list2):
    gene_pairs = []
    for elem1 in list1:
        for elem2 in list2:
            gene_pairs.append((elem1, elem2))
    return gene_pairs


def get_ortholog_pairs_by_species(t, sp1, sp2):
    def get_sps_info(t):
        sp_count = defaultdict(list)
        for i in t.get_leaf_names():
            sps = i.split("_")[0]
            sp_count[sps].append(i)
        return sp_count

    child1, child2 = t.get_children()
    dic1 = get_sps_info(child1)
    dic2 = get_sps_info(child2)

    a_sp1_gene_list = list(dic1[sp1])
    a_sp2_gene_list = list(dic1[sp2])
    b_sp1_gene_list = list(dic2[sp1])
    b_sp2_gene_list = list(dic2[sp2])

    gene_pairs = generate_combinations(a_sp1_gene_list, b_sp2_gene_list) + generate_combinations(b_sp1_gene_list, a_sp2_gene_list)
    return gene_pairs


def collect_speciation_pairs(
    Phylo_t: object,
    sp1,
    sp2,
    new_named_gene2gene_dic,
    processed_lines,
    written_results,
    pair_support=50,
) -> None:
    """Collect orthology-support pairs from speciation events.

    This preserves the legacy HaploFinder behavior for red labels while
    duplication-node detection itself is delegated to ``__init__.find_dup_node``.
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
        orthology = get_ortholog_pairs_by_species(spec_node, sp1, sp2)
        for gene1, gene2 in orthology:
            pair_node = Phylo_t.get_common_ancestor([gene1, gene2])
            if judge_support(pair_node.support, pair_support):
                gene_a = new_named_gene2gene_dic[gene1]
                gene_b = new_named_gene2gene_dic[gene2]
                result = f"{gene_a}\t{gene_b}\tred\n"
                token = f"{gene2}-r"
                if token not in written_results:
                    processed_lines.append(result)
                    written_results.add(token)


def find_independent_dup_nodes(dup_node_list) -> list:
    dup_node_list.sort(key=lambda x: len(x.get_leaf_names()))
    node_dic = {}

    for node in dup_node_list:
        if len(get_species_set(node)) == 1:
            continue
        else:
            tips = len(node)
            if tips in node_dic:
                node_dic[tips].append(node)
            else:
                node_dic[tips] = [node]
    if node_dic:
        min_key_value = min(node_dic.items(), key=lambda x: x[0])

        return min_key_value[1]
    else:
        return []


# ======================================================
# Section 6: Gene Conversion Detection
# ======================================================


def find_gene_conversion(dict_gd, dict_gff1, dict_gff2, lens_1, lens_2):
    chrs_combinations = [(chr_a, chr_b) for chr_a in [i[0] for i in lens_1]
                         for chr_b in [i[0] for i in lens_2]]

    block_list = []

    for chr_a, chr_b in chrs_combinations:
        digits_a = "".join(filter(str.isdigit, chr_a))
        digits_b = "".join(filter(str.isdigit, chr_b))
        if not digits_a or not digits_b:
            continue
        chr1 = int(digits_a)
        chr2 = int(digits_b)

        if chr1 * 2 == chr2:
            red_count = 0
            blue_count = 0

            for pair, color in dict_gd.items():
                gene_a, gene_b = pair.split(":")
                chr_gene_a = dict_gff1.get(gene_a, [None])[0]
                chr_gene_b = dict_gff2.get(gene_b, [None])[0]

                if chr_gene_a == chr_a and chr_gene_b == chr_b:
                    if color == "red":
                        red_count += 1
                    elif color == "blue":
                        blue_count += 1

            total = red_count + blue_count
            red_ratio = red_count / total if total > 0 else 0
            blue_ratio = blue_count / total if total > 0 else 0

            block_list.append(
                {
                    "chr_a": chr_a,
                    "chr_b": chr_b,
                    "red_count": red_count,
                    "blue_count": blue_count,
                    "red_ratio": red_ratio,
                    "blue_ratio": blue_ratio,
                }
            )
    filter_block_list = [
        block for block in block_list
        if block["blue_ratio"] > 0
    ]

    return filter_block_list


def find_gene_pair_info(gene_conversion_list, dict_gd, dict_gff1, dict_gff2, gd_pairs):
    sort_lst = []
    with open(f"gene_conversion_{gd_pairs}.txt", "w") as f:
        for block in gene_conversion_list:
            chr_a = block["chr_a"]
            chr_b = block["chr_b"]

            for pair, color in dict_gd.items():
                gene_a, gene_b = pair.split(":")
                chr_gene_a = dict_gff1.get(gene_a, [None])[0]
                chr_gene_b_info = dict_gff2.get(gene_b, [None, None, None])
                chr_gene_b = chr_gene_b_info[0]
                start_pos_b = chr_gene_b_info[1]

                if chr_gene_a == chr_a and chr_gene_b == chr_b and start_pos_b is not None:
                    start_pos_b_int = int(start_pos_b)
                    sort_lst.append((gene_a, gene_b, color, chr_gene_b, start_pos_b_int))

    sort_lst.sort(key=lambda x: (x[3], x[4]))

    new_lst = []
    with open(f"gene_conversion_{gd_pairs}.txt", "w") as f:
        for gene_a, gene_b, color, _, _ in sort_lst:
            f.write(f"{gene_a}\t{gene_b}\t{color}\n")
            new_lst.append((gene_a, gene_b, color))
    return new_lst


def find_conversion_zones_with_ids_to_file(data, dict_gff1, dict_gff2, output_file="gene_conversion.txt"):
    conversion_zones = []
    n = len(data)
    zone_id = 1
    i = 0

    while i < n - 1:
        if data[i][2] == "blue":
            start_red1 = i

            gene1, gene2, _ = data[start_red1]
            if gene1 not in dict_gff1 or gene2 not in dict_gff2:
                i += 1
                continue

            chrom1, chrom2 = dict_gff1[gene1][0], dict_gff2[gene2][0]

            while i < n and data[i][2] == "blue" and dict_gff1[data[i][0]][0] == chrom1 and dict_gff2[data[i][1]][0] == chrom2:
                i += 1
            end_red1 = i - 1

            if i < n and data[i][2] == "red" and dict_gff1[data[i][0]][0] == chrom1 and dict_gff2[data[i][1]][0] == chrom2:
                start_blue = i

                while i < n and data[i][2] == "red" and dict_gff1[data[i][0]][0] == chrom1 and dict_gff2[data[i][1]][0] == chrom2:
                    i += 1
                end_blue = i - 1

                if i < n and data[i][2] == "blue" and dict_gff1[data[i][0]][0] == chrom1 and dict_gff2[data[i][1]][0] == chrom2:
                    start_red2 = i

                    while i < n and data[i][2] == "blue" and dict_gff1[data[i][0]][0] == chrom1 and dict_gff2[data[i][1]][0] == chrom2:
                        i += 1
                    end_red2 = i - 1

                    conversion_zones.append((zone_id, start_red1, end_red2))
                    zone_id += 1

                    i = start_red2
                else:
                    i += 1
            else:
                i += 1
        else:
            i += 1

    with open(output_file, "w") as f:
        for zone in conversion_zones:
            zone_id, start, end = zone
            first_gene1, first_gene2, _ = data[start]
            contig1 = dict_gff1[first_gene1][0]
            contig2 = dict_gff2[first_gene2][0]
            f.write(f"# Conversion Zone {zone_id}: {contig1}&{contig2}\n")
            for j in range(start, end + 1):
                if j >= n:
                    break
                gene1, gene2, color = data[j]
                if gene1 in dict_gff1 and gene2 in dict_gff2:
                    chrom1, pos_start1, pos_end1 = dict_gff1[gene1][0], dict_gff1[gene1][1], dict_gff1[gene1][2]
                    chrom2, pos_start2, pos_end2 = dict_gff2[gene2][0], dict_gff2[gene2][1], dict_gff2[gene2][2]
                    position = f"{gene1}\t{pos_start1}\t{pos_end1}\t{gene2}\t{pos_start2}\t{pos_end2}"
                    f.write(f"{position}\t{color}\n")


# ======================================================
# Section 7: GD Pair Processing
# ======================================================


def process_gd_result(gf, imap, input_sps_tree, sp1, sp2, support, pair_support):
    tre_dic = read_and_return_dict(gf)
    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, taxa2voucher_dic = gene_id_transfer(imap)
    rename_sp1 = taxa2voucher_dic[sp1]
    rename_sp2 = taxa2voucher_dic[sp2]
    renamed_sptree = rename_input_tre(read_phylo_tree(input_sps_tree), taxa2voucher_dic)
    processed_lines = []
    written_results = set()
    for tre_ID, tre_path in tre_dic.items():
        Phylo_t0 = read_phylo_tree(tre_path)
        Phylo_t0.resolve_polytomy(recursive=True)
        Phylo_t0.sort_descendants()
        Phylo_t1 = rename_input_tre(Phylo_t0, gene2new_named_gene_dic)

        if len(get_species_set(Phylo_t1)) == 1:
            continue
        sps_tol = get_species_set(Phylo_t1)

        if len(sps_tol) == 2 and not {sp1, sp2}.issubset(sps_tol):
            continue

        annotate_gene_tree(Phylo_t1, renamed_sptree)
        dup_node_list_all = find_dup_node(
            Phylo_t1,
            renamed_sptree,
            gd_support=support,
            clade_support=support,
            dup_species_num=1,
            dup_species_percent=0,
            max_topology_distance=10**9,
        )
        dup_node_list = [
            node for node in dup_node_list_all
            if rename_sp1 in get_species_set(node) and rename_sp2 in get_species_set(node)
        ]
        collect_speciation_pairs(
            Phylo_t1,
            rename_sp1,
            rename_sp2,
            new_named_gene2gene_dic,
            processed_lines,
            written_results,
            pair_support,
        )

        for i in dup_node_list:
            if len(i) == 3:
                gene_pairs = get_ortholog_pairs_by_species(i, rename_sp1, rename_sp2)
                gd_clade1, gd_clade2 = i.get_children()
                gd_tips1 = set(gd_clade1.get_leaf_names())
                gd_tips2 = set(gd_clade2.get_leaf_names())

                for item in gene_pairs:
                    gene1, gene2 = item
                    item_set = set(item)
                    gene_a = new_named_gene2gene_dic[gene1]
                    gene_b = new_named_gene2gene_dic[gene2]

                    if item_set <= gd_tips1 or item_set <= gd_tips2:
                        result = f"{gene_a}\t{gene_b}\tred\n"
                        s = f"{gene2}-r"

                    else:
                        result = f"{gene_a}\t{gene_b}\tblue\n"
                        s = f"{gene2}-b"

                    if s not in written_results:
                        processed_lines.append(result)
                        written_results.add(s)

            elif len(i) == 4:
                gd_clade1, gd_clade2 = i.get_children()
                tips1 = set(get_species_set(gd_clade1))
                tips2 = set(get_species_set(gd_clade2))
                if len(tips1) == len(tips2) == 2:
                    gene_pairs = get_ortholog_pairs_by_species(i, rename_sp1, rename_sp2)
                    gd_tips1 = set(gd_clade1.get_leaf_names())
                    gd_tips2 = set(gd_clade2.get_leaf_names())
                    for item in gene_pairs:
                        gene1, gene2 = item
                        item_set = set(item)
                        gene_a = new_named_gene2gene_dic[gene1]
                        gene_b = new_named_gene2gene_dic[gene2]

                        if item_set <= gd_tips1 or item_set <= gd_tips2:
                            result = f"{gene_a}\t{gene_b}\tred\n"
                            s = f"{gene2}-r"

                        else:
                            result = f"{gene_a}\t{gene_b}\tblue\n"
                            s = f"{gene2}-b"

                        if s not in written_results:
                            processed_lines.append(result)
                            written_results.add(s)

    sorted_lines = sorted(processed_lines, key=lambda x: x.split("\t")[0])
    with open("color_label.txt", "w") as file:
        for line in sorted_lines:
            file.write(f"{line}")
    return sorted_lines


# ======================================================
# Section 8: FASTA Utilities
# ======================================================


def read_fasta(file_path):
    sequences = {}
    with open(file_path, "r") as f:
        current_id = None
        current_seq = []
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


def write_fasta(sequences, file_path):
    with open(file_path, "w") as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")


# ======================================================
# Section 9: Single-Copy Tree Extraction
# ======================================================


def get_single_copy_tree(tree, tree_id, tree_path, gene2new_named_gene_dic, renamed_length_dict=None):
    renamed_len_dic = {k: 0 for k in tree.get_leaf_names()}
    principal_gene_set, filtered_offcut_ev_seqs = offcut_tre(tree, renamed_len_dic)
    minor_orthologs = []
    minor_orthologs = iterator(
        filtered_offcut_ev_seqs,
        tree,
        gene2new_named_gene_dic,
        minor_orthologs,
        tree_path,
        renamed_length_dict,
    )
    ordered_name_OG_list = rename_OGs_tre_name(principal_gene_set, minor_orthologs, tree_id)
    sctree = []
    for tre_name, OG_set in ordered_name_OG_list:
        phylo_tree_0 = read_phylo_tree(tree_path)
        phylo_tree = root_tre_with_midpoint_outgroup(phylo_tree_0)
        rename_t = rename_input_tre(phylo_tree, gene2new_named_gene_dic)
        phylo_tree_OG_list = extract_tree(OG_set, rename_t)
        sctree.append(phylo_tree_OG_list)

    return sctree


# ======================================================
# Section 10: Hybrid Subgenome Assignment
# ======================================================


def assign_hybrid_subgenome(tree, hybrid_prefix, diploid_tags):
    distance_threshold = len(diploid_tags)

    result = {}
    subgenome_labels = [chr(ord("A") + i) for i in range(len(diploid_tags))]

    def is_hybrid(name):
        return hybrid_prefix in name

    def get_diploid_type(name):
        for tag in diploid_tags:
            if tag in name:
                return tag
        return None

    def calculate_distance_and_check(tree, distance_threshold, hybrid_leaves, diploid_leaves):
        if not hybrid_leaves or not diploid_leaves:
            return False

        hybrid_node = tree.get_common_ancestor(hybrid_leaves)
        diploid_node = tree.get_common_ancestor(diploid_leaves)

        if hybrid_node and diploid_node:
            node_distance = tree.get_distance(hybrid_node, diploid_node, topology_only=True)
            return node_distance <= distance_threshold

    def assign_subgenome(hybrid_list, diploid_types):
        if len(diploid_types) != 1:
            return False

        diploid_type = list(diploid_types)[0]
        if diploid_type not in diploid_tags:
            return False

        idx = diploid_tags.index(diploid_type)
        for h in hybrid_list:
            if is_hybrid(h):
                result[h] = subgenome_labels[idx]
        return True

    def traverse(node):
        if node.is_leaf():
            return
        if len(node.children) != 2:
            for child in node.children:
                traverse(child)
            return

        leaves1 = node.children[0].get_leaf_names()
        leaves2 = node.children[1].get_leaf_names()
        hyb1 = [x for x in leaves1 if is_hybrid(x)]
        hyb2 = [x for x in leaves2 if is_hybrid(x)]
        dip1 = [x for x in leaves1 if get_diploid_type(x)]
        dip2 = [x for x in leaves2 if get_diploid_type(x)]

        dip1_types = set(get_diploid_type(x) for x in dip1 if get_diploid_type(x))
        dip2_types = set(get_diploid_type(x) for x in dip2 if get_diploid_type(x))

        assigned = False

        if hyb1 and not dip1 and not hyb2 and dip2 and len(dip2_types) == 1:
            assigned = assign_subgenome(hyb1, dip2_types)

        elif hyb1 and not dip1 and hyb2 and dip2 and len(dip2_types) == 1:
            assigned = assign_subgenome(hyb1 + hyb2, dip2_types)

        elif hyb2 and not dip2 and not hyb1 and dip1 and len(dip1_types) == 1:
            assigned = assign_subgenome(hyb2, dip1_types)
        elif hyb2 and not dip2 and hyb1 and dip1 and len(dip1_types) == 1:
            assigned = assign_subgenome(hyb2 + hyb1, dip1_types)

        if not assigned:
            traverse(node.children[0])
            traverse(node.children[1])

    traverse(tree)

    for leaf in tree.get_leaf_names():
        if is_hybrid(leaf) and leaf not in result:
            result[leaf] = "unknown"

    return result


def split_sequences(input_GF_list, input_imap, hyb_sps, parental_sps, gff):
    tre_dic = read_and_return_dict(input_GF_list)
    gff_1, dict_gff1 = read_gff(gff)
    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, taxa2voucher_dic = gene_id_transfer(input_imap)
    hyb_sps_renamed = taxa2voucher_dic.get(hyb_sps, hyb_sps)
    parental_sps_renamed = [taxa2voucher_dic.get(s, s) for s in parental_sps]
    all_gene_labels = {}
    c = 0
    for tre_ID, tre_path in tre_dic.items():
        print(tre_ID)
        Phylo_t0 = read_phylo_tree(tre_path)
        Phylo_t0.resolve_polytomy(recursive=True)
        Phylo_t0.sort_descendants()
        Phylo_t1 = rename_input_tre(Phylo_t0, gene2new_named_gene_dic)
        leaf_names = Phylo_t1.get_leaf_names()
        leaf_species = set([g.split("_")[0] for g in leaf_names])
        if len(leaf_species) == 1:
            continue
        if not any(hyb_sps_renamed in leaf_name for leaf_name in leaf_names):
            continue
        a = assign_hybrid_subgenome(Phylo_t1, hyb_sps_renamed, parental_sps_renamed)
        b = {new_named_gene2gene_dic[k]: v for k, v in a.items()}

        for k1, v1 in b.items():
            if k1 in dict_gff1:
                c = dict_gff1[k1][0]
                is_correct_assignment = get_chromosome_subgenome(c)
                d = "unknown" if v1 == "unknown" else (v1 == is_correct_assignment)
                print(k1, v1, c, is_correct_assignment, d)

        print("-" * 30)


# ======================================================
# Section 11: Chromosome Mapping
# ======================================================


def get_chromosome_subgenome(chr_name):
    """
    Map chromosome identifiers to subgenome labels.

    Parameters
    ----------
    chr_name : str
        Chromosome identifier.

    Returns
    -------
    str
        Subgenome label or None if not assigned.

    Assumptions
    -----------
    Chromosome numbers encode subgenome A (1-10) and B (11-20).
    """
    match = re.search(r"(\d+)", chr_name)
    if match:
        chr_num = int(match.group(1))
        if 1 <= chr_num <= 10:
            return "A"
        elif 11 <= chr_num <= 20:
            return "B"
    return None


# ======================================================
# Section 12: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    gff1 = sys.argv[1]
    gff2 = sys.argv[2]
    lens1 = sys.argv[3]
    lens2 = sys.argv[4]
    blastp_pairs = sys.argv[5]
    synteny_pairs = sys.argv[6]
    spe1 = sys.argv[7]
    spe2 = sys.argv[8]
    num = int(sys.argv[9])
    gf = sys.argv[10]
    imap = sys.argv[11]
    target_chr1 = sys.argv[12]
    target_chr2 = sys.argv[13]
    size = float(sys.argv[14])

    # TODO: user will supplement process_blastp_result and parse_synteny_file
    process_blastp_pairs = process_blastp_result(blastp_pairs, num)
    alignments, alignment_scores = parse_synteny_file(synteny_pairs)
    process_synteny_pairs = assign_colors_by_alignment(alignments, alignment_scores)
    input_sps_tree = sys.argv[15]
    process_gd_pairs = process_gd_result(gf, imap, input_sps_tree, spe1, spe2, 50, 50)
    generate_dotplot(gff1, gff2, lens1, lens2, process_blastp_pairs, spe1, spe2, "blastp_pairs", target_chr1, target_chr2, size)
    print("-" * 30)
    generate_dotplot(gff1, gff2, lens1, lens2, process_synteny_pairs, spe1, spe2, "synteny_pairs", target_chr1, target_chr2, size)
    print("-" * 30)
    generate_dotplot(gff1, gff2, lens1, lens2, process_gd_pairs, spe1, spe2, "gd_pairs", target_chr1, target_chr2, size)
    total_pairs = [process_blastp_pairs, process_synteny_pairs, process_gd_pairs]
    # TODO: user will supplement process_total_color_list
    total_lst = process_total_color_list(total_pairs)
    generate_dotplot(gff1, gff2, lens1, lens2, total_lst, spe1, spe2, "total_pairs", target_chr1, target_chr2, size)
