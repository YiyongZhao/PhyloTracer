"""
Hybridization visualization utilities for the PhyloTracer pipeline.

This module parses HyDe outputs, computes summary statistics, and produces
annotated tree and heatmap figures for hybridization inference.
"""

import os
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from ete3 import NodeStyle, TextFace, Tree, TreeStyle
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from PIL import Image, ImageDraw, ImageFont

from phylotracer import (
    num_tre_node,
    realign_branch_length,
    rejust_root_dist,
)

# ======================================================
# Section 1: HyDe Output Parsing
# ======================================================


def parse_hyde_out(filename):
    """
    Parse HyDe output and return filtered result tuples.

    Parameters
    ----------
    filename : str
        Path to a HyDe output file.

    Returns
    -------
    list
        List of filtered HyDe entries.

    Assumptions
    -----------
    HyDe outputs are tab-delimited with a header beginning with 'P1'.
    """
    lst = []
    with open(filename, "r") as f:
        for i in f:
            i1 = i.strip().split("\t")
            if i1[0] == "P1":
                continue
            elif i1[5] == "nan":
                continue
            elif 0 < float(i1[5]) < 1:
                tup = [i1[0], i1[1], i1[2], i1[4], i1[5]]
                lst.append(tup)
    return lst


def calculate_three_tup(hyde_out_lst):
    """
    Group HyDe outputs by triplet key.

    Parameters
    ----------
    hyde_out_lst : list
        List of HyDe output tuples.

    Returns
    -------
    dict
        Mapping from triplet keys to lists of entries.

    Assumptions
    -----------
    Each entry has at least three taxon identifiers.
    """
    count_dic = defaultdict(list)
    for i in hyde_out_lst:
        key = "-".join(i[:3])
        count_dic[key].append(i)
    return count_dic


def get_hybrid_dic(summary_dic):
    """
    Group triplet keys by hybrid species.

    Parameters
    ----------
    summary_dic : dict
        Mapping from triplet keys to lists of entries.

    Returns
    -------
    dict
        Mapping from hybrid taxa to triplet key lists.

    Assumptions
    -----------
    Triplet keys are formatted as P1-Hybrid-P2.
    """
    hybrid_dic = {}
    for k, v in summary_dic.items():
        hybrid_tup = k.split("-")
        if hybrid_tup[1] in hybrid_dic:
            hybrid_dic[hybrid_tup[1]].append(k)
        else:
            hybrid_dic[hybrid_tup[1]] = [k]
    return hybrid_dic


def calculate_gamma(lst):
    """
    Compute mean gamma across a list of HyDe entries.

    Parameters
    ----------
    lst : list
        List of HyDe tuples.

    Returns
    -------
    float
        Mean gamma value.

    Assumptions
    -----------
    Gamma is stored at index 4 and is numeric when valid.
    """
    gamma = 0
    valid_count = 0
    for i in lst:
        if i[4] == "nan":
            continue
        elif i[4] == "-inf":
            continue
        else:
            num = float(i[4])
            gamma += num
            valid_count += 1
    return gamma / valid_count if valid_count else 0.0


def calculate_pvalue(lst):
    """
    Compute mean p-value across a list of HyDe entries.

    Parameters
    ----------
    lst : list
        List of HyDe tuples.

    Returns
    -------
    float
        Mean p-value.

    Assumptions
    -----------
    P-values are stored at index 3 and are numeric when valid.
    """
    pvalue = 0
    valid_count = 0
    for i in lst:
        if i[3] == "nan":
            continue
        elif i[3] == "-inf":
            continue
        else:
            num = float(i[3])
            pvalue += num
            valid_count += 1
    return pvalue / valid_count if valid_count else 0.0


# ======================================================
# Section 2: Tree Rendering Helpers
# ======================================================


def generate_tree_leaf(t, hybrid_sps, filename):
    """
    Render a tree with a highlighted hybrid leaf.

    Parameters
    ----------
    t : object
        Tree object to render.
    hybrid_sps : object
        Species name to highlight.
    filename : str
        Output file prefix.

    Returns
    -------
    None

    Assumptions
    -----------
    Tree nodes support ETE styling and faces.
    """
    nstyle = NodeStyle()
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1
    nstyle["vt_line_type"] = 0
    nstyle["hz_line_type"] = 0
    nstyle["size"] = 0
    nstyle["shape"] = "circle"
    nstyle["fgcolor"] = "black"

    leaf_names = t.get_leaf_names()
    max_l = max(len(name) for name in leaf_names)

    for i in t.traverse():
        i.set_style(nstyle)
        if i.is_leaf():
            leaf_name = i.name.ljust(max_l, "·")
            fgcolor = "red" if i.name == hybrid_sps else "black"
            i.add_face(
                TextFace(
                    leaf_name + " ",
                    fgcolor=fgcolor,
                    ftype="Arial",
                    fstyle="italic",
                ),
                column=1,
                position="aligned",
            )

    ts = TreeStyle()
    ts.scale = 10
    ts.show_leaf_name = False
    ts.show_scale = False
    realign_branch_length(t)
    rejust_root_dist(t)

    t.render(file_name=filename + "_img_faces.png", h=3200, tree_style=ts)


def generate_tree_node(t, node, filename):
    """
    Render a tree highlighting a target internal node and its clade.

    Parameters
    ----------
    t : object
        Tree object to render.
    node : object
        Target node to highlight.
    filename : str
        Output file prefix.

    Returns
    -------
    None

    Assumptions
    -----------
    Tree nodes support ETE styling and faces.
    """
    nodes = [i.name for i in node.traverse()]
    max_l = max([len(j) for j in t.get_leaf_names()])

    def set_node_style(n, color):
        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0
        nstyle["hz_line_type"] = 0
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["vt_line_color"] = color
        nstyle["hz_line_color"] = color
        n.set_style(nstyle)

    for i in t.traverse():
        if i.name == node.name:
            i.add_face(TextFace(i.name, fsize=4, fgcolor="blue"), column=1, position="branch-top")
        if i.name in nodes:
            set_node_style(i, "red")
        else:
            set_node_style(i, "black")

    for leaf in t:
        if leaf.name in node.get_leaf_names():
            leaf.add_face(
                TextFace(
                    leaf.name.ljust(max_l, "·"),
                    fgcolor="red",
                    fsize=10,
                    ftype="Arial",
                    fstyle="italic",
                ),
                column=1,
                position="aligned",
            )
        else:
            leaf.add_face(
                TextFace(
                    leaf.name.ljust(max_l, "·"),
                    fgcolor="black",
                    fsize=10,
                    ftype="Arial",
                    fstyle="italic",
                ),
                column=1,
                position="aligned",
            )

    ts = TreeStyle()
    ts.scale = 10
    ts.show_leaf_name = False

    ts.show_scale = False

    realign_branch_length(t)
    rejust_root_dist(t)

    t.render(file_name=filename + "_img_faces.png", h=3200, tree_style=ts)


# ======================================================
# Section 3: Summary-to-Matrix Utilities
# ======================================================


def from_summary_get_hyb_to_date(summary_dic, a_hyb_to_three_tup_list, leafs):
    """
    Convert summary dictionaries to matrices for visualization.

    Parameters
    ----------
    summary_dic : dict
        Mapping from triplet keys to HyDe output lists.
    a_hyb_to_three_tup_list : list
        List of triplet keys for a hybrid taxon.
    leafs : list
        Leaf names defining matrix order.

    Returns
    -------
    tuple
        (tup_df, gamma_df, pvalue_df) matrices.

    Assumptions
    -----------
    Leaf identifiers match triplet keys.
    """
    tup_df = pd.DataFrame(index=leafs, columns=leafs, data=0, dtype=float)
    gamma_df = pd.DataFrame(index=leafs, columns=leafs, data=0, dtype=float)
    pvalue_df = pd.DataFrame(index=leafs, columns=leafs, data=0, dtype=float)
    for i in a_hyb_to_three_tup_list:
        tup = i.split("-")
        p1, hyb, p2 = tup
        same_tups = summary_dic[i]
        num = len(same_tups)
        gamma = calculate_gamma(same_tups)
        pvalue = calculate_pvalue(same_tups)

        gamma_3 = round(gamma, 3)
        tup_df.loc[p1, p2] = num
        gamma_df.loc[p1, p2] = gamma_3
        pvalue_df.loc[p1, p2] = pvalue

    return tup_df, gamma_df, pvalue_df


# ======================================================
# Section 4: Colormap Utilities
# ======================================================


def hyde_visual_cmap():
    """
    Construct a custom colormap for HyDe visualizations.

    Returns
    -------
    ListedColormap
        Custom blue-to-red colormap.

    Assumptions
    -----------
    Colormap resolution is fixed at 10,000 levels.
    """
    cmap = []
    color = 1
    lucency = 0
    for each in range(5000):
        cmap.insert(each, np.array([0, 0, color, lucency]))
        color = color - (1 / 5000)
        lucency = lucency + (1 / 5000)
    color = 0
    lucency = 1
    for each in range(5000, 10000):
        cmap.insert(each, np.array([color, 0, 0, lucency]))
        color = color + (1 / 5000)
        lucency = lucency - (1 / 5000)
    newcmp = ListedColormap(cmap)
    return newcmp


def get_necmap():
    """
    Generate a combined lightened colormap for heatmaps.

    Returns
    -------
    LinearSegmentedColormap
        Combined colormap for heatmap rendering.

    Assumptions
    -----------
    Two base RGB colors define the gradient extremes.
    """

    def create_lightened_cmap(rgb_color, reverse=False):
        colors = [rgb_color]
        for i in range(50):
            lightened_color = tuple(
                (np.array(rgb_color) + (1 - np.array(rgb_color)) * (i / 50)).tolist()
            )
            colors.append(lightened_color)
        cmap_values = np.linspace(0, 1, len(colors))
        if reverse:
            colors = colors[::-1]
        lightened_cmap = LinearSegmentedColormap.from_list(
            "lightened_cmap",
            list(zip(cmap_values, colors)),
        )
        return lightened_cmap

    rgb_color1 = (33 / 255, 205 / 255, 67 / 255)
    rgb_color2 = (33 / 255, 171 / 255, 205 / 255)

    lightened_cmap1 = create_lightened_cmap(rgb_color1, reverse=True)
    lightened_cmap2 = create_lightened_cmap(rgb_color2)

    combined_cmap = LinearSegmentedColormap.from_list(
        "combined_cmap",
        np.vstack(
            (
                lightened_cmap1(np.linspace(0, 1, 256)),
                lightened_cmap2(np.linspace(0, 1, 256)),
            )
        ),
    )

    return combined_cmap


# ======================================================
# Section 5: Heatmap Construction
# ======================================================


def create_hot_map_node(summary_dic, hyb_dic, node, t, filename):
    """
    Create a heatmap for hybridization summaries at an internal node.

    Parameters
    ----------
    summary_dic : dict
        Mapping from triplet keys to HyDe output lists.
    hyb_dic : dict
        Mapping from hybrid taxa to triplet key lists.
    node : object
        Internal node defining the clade of interest.
    t : object
        Tree object for context.
    filename : str
        Output file prefix.

    Returns
    -------
    None

    Assumptions
    -----------
    Heatmap matrices are square and indexed by leaf names.
    """

    def average_dataframes(summary_tup_df):
        tup_result_df = pd.DataFrame(
            0,
            index=summary_tup_df[0].index,
            columns=summary_tup_df[0].columns,
        )
        for df in summary_tup_df:
            tup_result_df += df
        tup_result_df = tup_result_df.div(len(summary_tup_df))
        return tup_result_df

    node_s = node.get_leaf_names()
    leafs = t.get_leaf_names()
    df_lst = []
    for i in node_s:
        a_hyb_to_three_tup_list = hyb_dic[i]
        tup_df, gamma_df, filtered_gamma_df = from_summary_get_hyb_to_date(
            summary_dic,
            a_hyb_to_three_tup_list,
            leafs,
        )
        df_lst.append((tup_df, gamma_df, filtered_gamma_df))

    summary_tup_df = [d[0] for d in df_lst]
    summary_gamma_df = [d[1] for d in df_lst]
    summary_filter_df = [d[2] for d in df_lst]

    tup_result_df = average_dataframes(summary_tup_df).astype(int)
    gamma_result_df = average_dataframes(summary_gamma_df)
    filter_result_df = average_dataframes(summary_filter_df)

    for leaf in node_s:
        tup_result_df.loc[:, leaf] = 0
        tup_result_df.loc[leaf] = 0
        gamma_result_df.loc[:, leaf] = 0
        gamma_result_df.loc[leaf] = 0

    mask_upper = np.triu(np.ones_like(gamma_result_df, dtype=bool), k=1)
    tup_result_df[mask_upper] = 0
    gamma_result_df[mask_upper] = 0

    tup_annot = tup_result_df.astype(str).where(tup_result_df != 0, other="")
    gamma_annot = gamma_result_df.map(lambda x: f"{x:.3f}" if x != 0 else "")

    fig = plt.figure(figsize=(30, 30))

    border_width = 0.001
    ax_size = [0 + border_width, 0 + border_width, 1 - 2 * border_width, 1 - 2 * border_width]
    ax = fig.add_axes(ax_size)

    newcmp = hyde_visual_cmap()

    sns.heatmap(
        tup_result_df,
        annot=tup_annot,
        fmt="",
        cmap="Greys",
        ax=ax,
        annot_kws={"color": "#FFFFFF", "size": 60, "va": "top"},
        xticklabels=False,
        yticklabels=False,
        cbar=False,
        linewidths=1.5,
        linecolor="black",
        square=True,
    )

    sns.heatmap(
        gamma_result_df,
        annot=gamma_annot,
        fmt="",
        cmap="Greys",
        ax=ax,
        annot_kws={"color": "#F5EF70", "size": 60, "va": "bottom"},
        xticklabels=False,
        yticklabels=False,
        cbar=False,
        linewidths=1.5,
        linecolor="black",
        square=True,
    )

    newcmp = hyde_visual_cmap()

    sns.heatmap(
        gamma_result_df,
        annot=False,
        cmap=newcmp,
        ax=ax,
        xticklabels=False,
        yticklabels=False,
        cbar=False,
        cbar_kws={"shrink": 0.6},
        vmin=0,
        vmax=1,
        linewidths=1.5,
        linecolor="black",
        square=True,
    )

    m = ax.imshow(gamma_result_df, norm=colors.Normalize(vmin=0, vmax=1), cmap=newcmp)
    position = fig.add_axes([0.9, 0.2, 0.05, 0.7])
    cbar = plt.colorbar(m, cax=position)
    cbar.ax.tick_params(labelsize=40)
    plt.savefig(filename + "_hotmap.png", dpi=200)
    plt.cla()
    plt.close("all")


def create_hot_map_leaf(summary_dic, a_hyb_to_three_tup_list, t, filename):
    """
    Create a heatmap for hybridization summaries at the leaf level.

    Parameters
    ----------
    summary_dic : dict
        Mapping from triplet keys to HyDe output lists.
    a_hyb_to_three_tup_list : list
        Triplet keys for a given hybrid.
    t : object
        Tree object defining leaf order.
    filename : str
        Output file prefix.

    Returns
    -------
    None

    Assumptions
    -----------
    Leaf names in the tree match those in triplet keys.
    """
    sp = t.get_leaf_names()
    tup_df = pd.DataFrame(index=sp, columns=sp, data=0, dtype=float)
    gamma_df = pd.DataFrame(index=sp, columns=sp, data=0, dtype=float)
    pvalue_df = pd.DataFrame(index=sp, columns=sp, data=0, dtype=float)

    for i in a_hyb_to_three_tup_list:
        tup = i.split("-")
        p1, hyb, p2 = tup
        if p1 not in sp or p2 not in sp:
            print(f"Warning: {p1} or {p2} is not in sp.")
            continue
        same_tups = summary_dic[i]
        num = len(same_tups)
        gamma = calculate_gamma(same_tups)
        pvalue = calculate_pvalue(same_tups)

        tup_df.at[p1, p2] = num
        gamma_df.at[p1, p2] = round(gamma, 3)
        pvalue_df.at[p1, p2] = pvalue

    fig = plt.figure(figsize=(30, 30))

    border_width = 0.001
    ax_size = [0 + border_width, 0 + border_width, 1 - 2 * border_width, 1 - 2 * border_width]
    ax = fig.add_axes(ax_size)

    processed_pairs = set()
    for p1 in sp:
        for p2 in sp:
            if p1 != p2:
                pair_key = tuple(sorted([p1, p2]))
                if pair_key in processed_pairs:
                    continue
                processed_pairs.add(pair_key)

                pvalue_p1_p2 = pvalue_df.loc[p1, p2]
                pvalue_p2_p1 = pvalue_df.loc[p2, p1]

                if pvalue_p1_p2 > pvalue_p2_p1:
                    gamma_df.loc[p2, p1] = 0
                else:
                    gamma_df.loc[p1, p2] = 0

    tup_annot = tup_df.astype(str).where(tup_df != 0, other="")
    gamma_annot = gamma_df.map(lambda x: f"{x:.3f}" if x != 0 else "")
    sns.heatmap(
        tup_df,
        annot=tup_annot,
        fmt="",
        cmap="Greys",
        ax=ax,
        annot_kws={"color": "#FFFFFF", "va": "top"},
        xticklabels=False,
        yticklabels=False,
        cbar=False,
        linewidths=1.5,
        linecolor="black",
        square=True,
    )

    sns.heatmap(
        gamma_df,
        annot=gamma_annot,
        fmt="",
        cmap="Greys",
        ax=ax,
        annot_kws={"color": "#F5EF70", "va": "bottom"},
        xticklabels=False,
        yticklabels=False,
        cbar=False,
        linewidths=1.5,
        linecolor="black",
        square=True,
    )

    newcmp = hyde_visual_cmap()

    sns.heatmap(
        gamma_df,
        annot=False,
        cmap=newcmp,
        ax=ax,
        xticklabels=False,
        yticklabels=False,
        cbar=False,
        cbar_kws={"shrink": 0.6},
        vmin=0,
        vmax=1,
        linewidths=1.5,
        linecolor="black",
        square=True,
    )

    m = ax.imshow(gamma_df, norm=colors.Normalize(vmin=0, vmax=1), cmap=newcmp)
    position = fig.add_axes([0.9, 0.2, 0.05, 0.7])
    cbar = plt.colorbar(m, cax=position)
    cbar.ax.tick_params(labelsize=40)
    cbar.set_label("Heatmap of y Values", fontsize=20, labelpad=2)

    plt.savefig(filename + "_hotmap.png", dpi=300)
    plt.cla()
    plt.close("all")


def combine_fig(hybrid_species):
    """
    Combine tree and heatmap images into a single figure.

    Parameters
    ----------
    hybrid_species : str
        Hybrid species identifier used for file naming.

    Returns
    -------
    None

    Assumptions
    -----------
    Required image files exist and can be loaded by PIL.
    """
    treepic = Image.open(f"{hybrid_species}_img_faces.png")
    treepic_size = treepic.size

    treepic_rotate = treepic.rotate(90, expand=1)
    rotate_size = treepic_rotate.size

    min_width = max(
        treepic_size[0] + treepic_size[1] + 60,
        treepic_size[0] + rotate_size[0] + 60,
    )
    min_height = max(
        treepic_size[1] + rotate_size[1] + 60,
        treepic_size[0] + treepic_size[1] + 60,
    )
    combine_fig_size = max(min_width, min_height)

    combine = Image.new("RGB", (combine_fig_size, combine_fig_size), "#FFFFFF")

    combine.paste(treepic, (40, 40))

    rotate_x = min(treepic_size[0] + 20, combine_fig_size - rotate_size[0] - 20)
    rotate_y = min(treepic_size[1] + 20, combine_fig_size - rotate_size[1] - 20)
    combine.paste(treepic_rotate, (rotate_x, rotate_y))

    hotpic = Image.open(f"{hybrid_species}_hotmap.png")
    hotpic.thumbnail((treepic_size[1], treepic_size[1]))
    hotmap_x = min(treepic_size[0] + 20, combine_fig_size - hotpic.size[0] - 20)
    hotmap_y = 40
    combine.paste(hotpic, (hotmap_x, hotmap_y))

    draw = ImageDraw.Draw(combine)

    base_font_size = max(16, int(combine_fig_size / 50))
    title_font_size = max(18, int(combine_fig_size / 45))
    legend_font_size = max(14, int(combine_fig_size / 55))

    try:
        font_base = ImageFont.truetype("/System/Library/Fonts/Arial.ttf", base_font_size)
        font_title = ImageFont.truetype("/System/Library/Fonts/Arial.ttf", title_font_size)
        font_legend = ImageFont.truetype("/System/Library/Fonts/Arial.ttf", legend_font_size)
    except Exception:
        font_base = ImageFont.load_default()
        font_title = ImageFont.load_default()
        font_legend = ImageFont.load_default()

    max_width = combine_fig_size
    max_height = combine_fig_size

    y_text = "y (hybridization)"
    y_text_x = max(10, min(50, max_width - 200))
    y_text_y = 15
    draw.text((y_text_x, y_text_y), y_text, fill="black", font=font_title)

    oney_text = "1-y (complement)"
    text_width_estimate = len(oney_text) * (title_font_size * 0.6)
    oney_text_x = max(10, min(rotate_x + 10, max_width - text_width_estimate - 20))
    oney_text_y = min(rotate_y + treepic_rotate.size[1] + 5, max_height - title_font_size - 10)

    legend_x = oney_text_x
    legend_spacing = max(18, int(legend_font_size * 1.3))
    legend_end_y = oney_text_y - 15
    legend_start_y = legend_end_y - (legend_spacing * 3)

    legend_text_width = max(
        len("Red: hybrid"),
        len("Yellow: y values"),
        len("White: combinations"),
    ) * (legend_font_size * 0.6)
    legend_x = max(10, min(legend_x, max_width - legend_text_width - 20))

    if legend_start_y > 0:
        draw.text((legend_x, legend_start_y), "Red: hybrid", fill="red", font=font_legend)
        draw.text(
            (legend_x, legend_start_y + legend_spacing),
            "Yellow: y values",
            fill="yellow",
            font=font_legend,
        )
        draw.text(
            (legend_x, legend_start_y + legend_spacing * 2),
            "White: number of hybridization combinations",
            fill="black",
            font=font_legend,
        )

    draw.text((oney_text_x, oney_text_y), oney_text, fill="black", font=font_title)

    combine.save(hybrid_species + ".png")


# ======================================================
# Section 6: Main Visualization Pipelines
# ======================================================


def hyde_visual_leaf_main(out_file_name, sptree):
    """
    Generate hybridization visualizations centered on leaves.

    Parameters
    ----------
    out_file_name : str
        Path to HyDe output file.
    sptree : object
        Species tree object.

    Returns
    -------
    None

    Assumptions
    -----------
    Output files are writable in the current working directory.
    """
    out1 = parse_hyde_out(out_file_name)
    result1 = calculate_three_tup(out1)
    hybrid_dic1 = get_hybrid_dic(result1)
    for k, v in hybrid_dic1.items():
        print(f"{k} is processing")
        t1 = sptree.copy()
        generate_tree_leaf(t1, k, k)
        create_hot_map_leaf(result1, v, t1, k)
        combine_fig(k)
        os.remove(f"{k}_hotmap.png")
        os.remove(f"{k}_img_faces.png")


def hyde_visual_node_main(out_file_name, sptree):
    """
    Generate hybridization visualizations centered on internal nodes.

    Parameters
    ----------
    out_file_name : str
        Path to HyDe output file.
    sptree : object
        Species tree object.

    Returns
    -------
    None

    Assumptions
    -----------
    Internal nodes are numbered prior to visualization.
    """
    out1 = parse_hyde_out(out_file_name)
    result1 = calculate_three_tup(out1)
    hybrid_dic1 = get_hybrid_dic(result1)

    num_tre_node(sptree)
    nodes = [i for i in sptree.traverse() if not i.is_leaf()]

    for node in nodes:
        if node.is_root():
            continue
        else:
            print(f"{node.name} is processing")
            t1 = sptree.copy()
            generate_tree_node(t1, node, node.name)
            create_hot_map_node(result1, hybrid_dic1, node, t1, node.name)
            combine_fig(node.name)
            os.remove(f"{node.name}_hotmap.png")
            os.remove(f"{node.name}_img_faces.png")


# ======================================================
# Section 7: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Visualize hybridization (HyDe) results on a species tree.",
    )
    parser.add_argument(
        "species_tree_file",
        help="Path to the species tree file (Newick format).",
    )
    parser.add_argument(
        "hyde_output_file",
        help="Path to the HyDe output file.",
    )
    parser.add_argument(
        "--mode",
        choices=["leaf", "node"],
        default="leaf",
        help="Visualization mode: 'leaf' for leaf-centered, 'node' for node-centered (default: leaf).",
    )
    args = parser.parse_args()

    sptree = Tree(args.species_tree_file)
    if args.mode == "leaf":
        hyde_visual_leaf_main(args.hyde_output_file, sptree)
    else:
        hyde_visual_node_main(args.hyde_output_file, sptree)
