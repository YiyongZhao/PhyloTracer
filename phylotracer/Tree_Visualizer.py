"""
Tree rendering and annotation utilities for the PhyloTracer pipeline.

This module provides visualization helpers for gene and species trees,
including branch-length realignment, duplication annotation, and
multi-category tip labeling for publication-quality figures.
"""
from __future__ import annotations

import os
import re
import shutil
import string

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ete3 import PhyloTree, Tree
try:
    from ete3 import NodeStyle, RectFace, TextFace, TreeStyle
except ImportError:
    NodeStyle = None
    RectFace = None
    TextFace = None
    TreeStyle = None
from tqdm import tqdm

from phylotracer import (
    find_dup_node,
    gene_id_transfer,
    is_rooted,
    num_tre_node,
    read_and_return_dict,
    read_phylo_tree,
    realign_branch_length,
    rejust_root_dist,
    rename_input_tre,
    root_tre_with_midpoint_outgroup,
)


# ======================================================
# Section 2: Duplication and TreeStyle Utilities
# ======================================================


def dup_nodeids_from_numbered_gfs(
    phylo_t: object,
    species_tree: object,
) -> tuple:
    """
    Number nodes and infer duplication events in a gene tree.

    Parameters
    ----------
    phylo_t : object
        ETE gene tree to be rooted and numbered.
    species_tree : object
        ETE species tree used for duplication inference.

    Returns
    -------
    tuple
        A tuple of (numbered_tree, duplication_node_list).

    Assumptions
    -----------
    Midpoint rooting is appropriate if the tree is not rooted.
    """
    if not is_rooted(phylo_t):
        phylo_t = root_tre_with_midpoint_outgroup(phylo_t)
    phylo_t = num_tre_node(phylo_t)
    dup_node_list = find_dup_node(phylo_t, species_tree)
    phylo_t1 = Tree(phylo_t.write())
    num_tre_node(phylo_t1)
    return phylo_t1, dup_node_list


def create_tree_style(tree_style: str, tree_id: str, visual: bool = False) -> object:
    """
    Create a TreeStyle configuration for gene tree rendering.

    Parameters
    ----------
    tree_style : str
        Display mode for the tree (e.g., "circular", "rectangular").
    tree_id : str
        Tree identifier displayed in the title.
    visual : bool, optional
        Whether to include duplication event legend entries.

    Returns
    -------
    object
        Configured TreeStyle object.

    Assumptions
    -----------
    TreeStyle and TextFace are available from the ETE toolkit.
    """
    ts = TreeStyle()
    ts.title.add_face(TextFace("★", fgcolor="white", ftype="Arial"), column=0)
    ts.title.add_face(TextFace(tree_id, ftype="Arial"), column=1)
    if visual:
        ts.title.add_face(TextFace("★", fgcolor="red", ftype="Arial"), column=0)
        ts.title.add_face(
            TextFace("Interspecific gene duplication event", ftype="Arial"),
            column=1,
        )
        ts.title.add_face(TextFace("★", fgcolor="blue", ftype="Arial"), column=0)
        ts.title.add_face(
            TextFace("Intraspecific gene duplication event", ftype="Arial"),
            column=1,
        )
    ts.mode = tree_style
    ts.scale = 20
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.extra_branch_line_type = 0
    ts.extra_branch_line_color = "black"
    ts.branch_vertical_margin = -1

    return ts


def set_node_style(
    node: object,
    dup_node_name_list: list,
    visual: bool = False,
) -> object:
    """
    Configure node styling and optionally mark duplication events.

    Parameters
    ----------
    node : object
        ETE node to be styled.
    dup_node_name_list : list
        List of node names that represent duplication events.
    visual : bool, optional
        Whether to add duplication event markers.

    Returns
    -------
    object
        Configured NodeStyle object.

    Assumptions
    -----------
    Species labels can be extracted by ``get_species_list``.
    """
    nstyle = NodeStyle()
    splist = set(get_species_list(node))
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1
    nstyle["vt_line_type"] = 0
    nstyle["hz_line_type"] = 0
    nstyle["size"] = 0
    nstyle["shape"] = "circle"
    nstyle["fgcolor"] = "black"

    if visual:
        if node.name in dup_node_name_list and len(splist) == 1:
            node.add_face(
                TextFace("★", fsize=7, fgcolor="blue", ftype="Arial"),
                column=1,
                position="branch-top",
            )
        elif node.name in dup_node_name_list and len(splist) != 1:
            node.add_face(
                TextFace("★", fsize=7, fgcolor="red", ftype="Arial"),
                column=1,
                position="branch-top",
            )

    node.set_style(nstyle)


# ======================================================
# Section 3: TreeStyle Assembly
# ======================================================


def get_treestyle(
    Phylo_t: object,
    tree_style: str,
    tre_ID: str,
    visual: bool = False,
    species_tree: object = None,
) -> object:
    """
    Build a styled tree and TreeStyle for rendering.

    Parameters
    ----------
    Phylo_t : object
        ETE gene tree object to be styled.
    tree_style : str
        Display mode for the tree.
    tre_ID : str
        Tree identifier for labeling.
    visual : bool, optional
        Whether to include duplication event markers.
    species_tree : object, optional
        Species tree used for duplication inference.

    Returns
    -------
    object
        A tuple of (styled_tree, tree_style_object).

    Assumptions
    -----------
    Duplication inference requires a valid species tree.
    """
    Phylo_t1, dup_node_list = dup_nodeids_from_numbered_gfs(
        Phylo_t,
        species_tree,
    )
    ts = create_tree_style(tree_style, tre_ID, visual)
    dup_node_name_list = [node.name for node in dup_node_list]
    for node in Phylo_t1.traverse():
        set_node_style(node, dup_node_name_list, visual)

    return Phylo_t1, ts


# ======================================================
# Section 4: Color Mapping Utilities
# ======================================================


def get_color_dict(dictionary: dict) -> dict:
    """
    Generate a color dictionary for categorical labels.

    Parameters
    ----------
    dictionary : dict
        Mapping from identifiers to category labels.

    Returns
    -------
    dict
        Mapping from identifiers to ``category@hexcolor`` strings.

    Assumptions
    -----------
    Categories are finite and can be mapped to a rainbow colormap.
    """
    colormap = plt.get_cmap("rainbow")
    unique_values = set(dictionary.values())
    colors_list = [
        colors.rgb2hex(colormap(i))
        for i in np.linspace(0, 1, len(unique_values))
    ]
    color_dict = dict(zip(unique_values, colors_list))
    sps_color_list = {
        k: f"{v}@{color_dict.get(v)}"
        for k, v in dictionary.items()
        if v in color_dict
    }
    return sps_color_list


def generate_color_dict(gene_categories: list[dict]) -> list[dict]:
    """
    Generate color dictionaries for multiple gene category mappings.

    Parameters
    ----------
    gene_categories : list[dict]
        List of dictionaries mapping gene identifiers to categories.

    Returns
    -------
    list[dict]
        List of color-mapped dictionaries corresponding to each input.

    Assumptions
    -----------
    Each category dictionary is compatible with ``get_color_dict``.
    """
    return [get_color_dict(category_dict) for category_dict in gene_categories]


def generate_string(index: int) -> str:
    """
    Generate a short suffix label for categorical legends.

    Parameters
    ----------
    index : int
        Index to convert into an alphanumeric suffix.

    Returns
    -------
    str
        Generated suffix string.

    Assumptions
    -----------
    Alphabetic labels are used up to 52 entries before reuse.
    """
    letters = list(string.ascii_uppercase) + list(string.ascii_lowercase)
    if index < len(letters):
        return "@" + letters[index]
    else:
        first_letter = letters[(index - len(letters)) // 52]
        second_letter = letters[(index - len(letters)) % 52]
        return "@" + first_letter + second_letter


def get_new_sorted_dict(gene2fam: dict) -> dict:
    """
    Create a sorted family-to-color mapping with unique suffixes.

    Parameters
    ----------
    gene2fam : dict
        Mapping from gene identifiers to family labels.

    Returns
    -------
    dict
        Mapping from family labels to ``hexcolor@suffix`` strings.

    Assumptions
    -----------
    Family labels are suitable for deterministic sorting.
    """
    uniq_fam = set(get_color_dict(gene2fam).values())
    fam2color = {i.split("@")[0]: i.split("@")[-1] for i in uniq_fam}
    sorted_dict = dict(sorted(fam2color.items(), key=lambda x: x[0], reverse=False))
    for index, key in enumerate(sorted_dict.keys()):
        suffix = generate_string(index)
        sorted_dict[key] += suffix
    return sorted_dict


def fuzzy_match(search_string: str, key: str) -> re.Match:
    """
    Perform a regex-based fuzzy match between a pattern and a key.

    Parameters
    ----------
    search_string : str
        Regex pattern used for matching.
    key : str
        Candidate string to match.

    Returns
    -------
    re.Match
        Match object if found, otherwise None.

    Assumptions
    -----------
    Regex patterns are valid for the Python ``re`` module.
    """
    return re.search(search_string, key)


# ======================================================
# Section 5: Tip Annotation and Heatmap Rendering
# ======================================================


def tips_mark(
    Phylo_t1: object,
    voucher2taxa_dic: dict,
    color_dicts: list,
    sps_color_dict: dict,
    tre_ID,
    ts,
    new_named_gene2gene_dic: dict,
    dir_path,
    gene_color_dict=None,
    df=None,
) -> object:
    """
    Annotate tree tips with categorical colors and optional heatmap values.

    Parameters
    ----------
    Phylo_t1 : object
        ETE tree object to annotate.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.
    color_dicts : list
        List of color dictionaries for gene categories.
    sps_color_dict : dict
        Color dictionary for species labels.
    tre_ID : object
        Tree identifier used in output filenames.
    ts : object
        TreeStyle configuration for rendering.
    new_named_gene2gene_dic : dict
        Mapping from renamed gene identifiers to original gene identifiers.
    dir_path : object
        Output directory for rendered PDFs.
    gene_color_dict : dict, optional
        Color dictionary for gene families.
    df : object, optional
        DataFrame containing heatmap values for tips.

    Returns
    -------
    object
        Rendered tree output from ``Tree.render``.

    Assumptions
    -----------
    The tree is compatible with ETE faces and rendering APIs.
    """

    def add_face_to_node(
        node: object,
        face: object,
        column: int,
        position: str = "aligned",
    ) -> None:
        """
        Add a face to a node only once at a given column and position.

        Parameters
        ----------
        node : object
            Tree node receiving the face.
        face : object
            Face object to be attached.
        column : int
            Column index for aligned faces.
        position : str, optional
            Face position relative to the node.

        Returns
        -------
        None

        Assumptions
        -----------
        Face identity can be tracked by node, column, and position.
        """
        if (node, column, position) not in faces_added:
            node.add_face(face, column=column, position=position)
            faces_added.add((node, column, position))

    def generate_face_mark(
        node: object,
        species: str,
        column: int,
        color_dict: dict,
    ) -> None:
        """
        Add a colored categorical face for a given species.

        Parameters
        ----------
        node : object
            Tree node to annotate.
        species : str
            Species label for the node.
        column : int
            Column index for the face.
        color_dict : dict
            Mapping from species to color-coded labels.

        Returns
        -------
        None

        Assumptions
        -----------
        Color dictionary values follow the ``label@hexcolor`` format.
        """
        if (node, column, "aligned") not in faces_added:
            if species in color_dict:
                color = color_dict[species].split("@")[ -1]
                face = TextFace(
                    "  ▐" + "  " + color_dict[species].split("@")[0],
                    fgcolor=color,
                    ftype="Arial",
                )
                add_face_to_node(node, face, column, position="aligned")
            else:
                color = "white"
                face = TextFace(
                    "  ▐" + "  " + "" * len(species),
                    fgcolor=color,
                    ftype="Arial",
                )
                add_face_to_node(node, face, column, position="aligned")

    def add_species_face(
        node: object,
        gene: str,
        species: str,
        sps_color_dict: dict,
    ) -> None:
        """
        Add gene and species labels with species-based coloring.

        Parameters
        ----------
        node : object
            Tree node to annotate.
        gene : str
            Gene label for the node.
        species : str
            Species label for the node.
        sps_color_dict : dict
            Species color dictionary.

        Returns
        -------
        None

        Assumptions
        -----------
        Species labels are present in ``sps_color_dict``.
        """
        if species in sps_color_dict:
            color = sps_color_dict[species].split("@")[ -1]
            gene_face = TextFace(
                " " + gene,
                fgcolor=color,
                ftype="Arial",
                fstyle="italic",
            )
            node.add_face(gene_face, column=-1)
        if species in sps_color_dict:
            color = sps_color_dict[species].split("@")[ -1]
            species_face = TextFace(
                "  ▐" + "  " + sps_color_dict[species].split("@")[0],
                fgcolor=color,
                ftype="Arial",
                fstyle="italic",
            )
            add_face_to_node(node, species_face, 0, position="aligned")

    def add_gene_face(node: object, gene: str, column: int) -> None:
        """
        Add a gene family face for matched gene patterns.

        Parameters
        ----------
        node : object
            Tree node to annotate.
        gene : str
            Gene label for matching.
        column : int
            Column index for the face.

        Returns
        -------
        None

        Assumptions
        -----------
        ``gene_color_dict`` uses ``label@hexcolor`` values.
        """
        matched_key = None
        matched_value = None
        for key in gene_color_dict:
            if fuzzy_match(key, gene):
                matched_key = key
                matched_value = gene_color_dict.get(key)
                break
        if matched_value:
            color = matched_value.split("@")[ -1]
            face5 = TextFace(
                "  ▐" + "  " + matched_value.split("@")[0],
                fgcolor=color,
                ftype="Arial",
            )
            add_face_to_node(node, face5, column, position="aligned")

    def get_color(value: float) -> str:
        """
        Map a numeric value to a discrete heatmap color.

        Parameters
        ----------
        value : float
            Numeric value for color mapping.

        Returns
        -------
        str
            Hex color code for the given value.

        Assumptions
        -----------
        Values outside defined bins use the highest intensity color.
        """
        if np.isnan(value):
            return "white"
        else:
            if 0 <= value <= 5:
                return "#006599"
            elif 5 < value <= 10:
                return "#408ca6"
            elif 10 < value <= 15:
                return "#7fb2b2"
            elif 15 < value <= 20:
                return "#bfd9bf"
            elif 20 < value <= 25:
                return "#ffffcc"
            elif 25 < value <= 30:
                return "#f7deab"
            elif 30 < value <= 35:
                return "#eebc88"
            elif 35 < value <= 40:
                return "#e69966"
            elif 40 < value <= 45:
                return "#dc7845"
            elif 45 < value <= 50:
                return "#d55623"
            else:
                return "#cc3300"

    def add_heat_map_to_node(
        tree: object,
        df: pd.DataFrame,
        new_named_gene2gene_dic: dict,
        start_col: int,
    ) -> None:
        """
        Add heatmap faces to tree tips based on a DataFrame.

        Parameters
        ----------
        tree : object
            Tree whose tips are annotated.
        df : pd.DataFrame
            DataFrame with rows keyed by gene labels.
        new_named_gene2gene_dic : dict
            Mapping from node names to gene names.
        start_col : int
            Starting column for heatmap faces.

        Returns
        -------
        None

        Assumptions
        -----------
        DataFrame indices can be matched to gene labels by regex.
        """
        columns = df.columns.tolist()
        for node in tree:
            gene = new_named_gene2gene_dic[node.name]
            matched_key = None
            for key in df.index:
                if fuzzy_match(key, gene):
                    matched_key = key
                    break
            for ind, col_name in enumerate(columns):
                col_idx = start_col + ind
                if node.is_leaf() and matched_key in df.index:
                    color = get_color(df[col_name][matched_key])
                else:
                    color = "#F5F5F5"
                face = RectFace(width=10, height=10, fgcolor=color, bgcolor=color)
                node.add_face(face, column=col_idx, position="aligned")

    def add_header_to_tree(ts: any, df: pd.DataFrame, new_start: int) -> None:
        """
        Add rotated column headers for heatmap values.

        Parameters
        ----------
        ts : any
            TreeStyle object receiving the header.
        df : pd.DataFrame
            DataFrame providing column labels.
        new_start : int
            Starting column for header placement.

        Returns
        -------
        None

        Assumptions
        -----------
        Column labels are short enough to be displayed vertically.
        """
        labels = df.columns.to_list()
        for ind, i in enumerate(labels):
            face = TextFace(" " + i, fgcolor="black", ftype="Arial", fsize=9)
            face.rotation = -90
            face.vt_align = 2
            ts.aligned_header.add_face(face, new_start + ind)

    def add_color_bar(ts: any) -> None:
        """
        Add a heatmap color bar legend to the tree visualization.

        Parameters
        ----------
        ts : any
            TreeStyle object receiving the legend.

        Returns
        -------
        None

        Assumptions
        -----------
        Heatmap bins are fixed to the predefined value ranges.
        """
        ts.legend.add_face(TextFace(" "), column=0)
        bar_face = TextFace("Color Bar ", ftype="Arial")
        ts.legend.add_face(bar_face, column=0)
        cols = [
            "#006599",
            "#408ca6",
            "#7fb2b2",
            "#bfd9bf",
            "#ffffcc",
            "#f7deab",
            "#eebc88",
            "#e69966",
            "#dc7845",
            "#d55623",
            "#cc3300",
        ]
        bounds = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        col_dic = dict(zip(bounds, cols))
        n = 1
        for k, v in col_dic.items():
            colorbar_face = RectFace(width=20, height=20, fgcolor=v, bgcolor=v)
            ts.legend.add_face(TextFace(" " + str(k)), column=n)
            ts.legend.add_face(colorbar_face, column=n)
            n += 1
        ts.legend_position = 2

    faces_added = set()

    for node in Phylo_t1.traverse():
        if node.is_leaf():
            gene = new_named_gene2gene_dic[node.name]
            species = voucher2taxa_dic[node.name.split("_")[0]]
            rename_species = node.name.split("_")[0]
            add_species_face(node, gene, rename_species, sps_color_dict)
            column = 1
            for color_dict in color_dicts:
                generate_face_mark(node, gene, column, color_dict)
                column += 1

            if gene_color_dict is not None and gene in gene_color_dict:
                add_gene_face(node, gene, column)
            else:
                face_placeholder = TextFace("  ▐", fgcolor="white", ftype="Arial")
                add_face_to_node(node, face_placeholder, column, position="aligned")
            column += 1
    if df is not None:
        add_heat_map_to_node(Phylo_t1, df, new_named_gene2gene_dic, column)
        add_header_to_tree(ts, df, column)
        add_color_bar(ts)

    return Phylo_t1.render(
        f'{dir_path}{tre_ID}.pdf',
        w=210,
        units="mm",
        tree_style=ts,
    )


# ======================================================
# Section 6: Gene Family and Species Tree Mapping
# ======================================================


def get_matched_value(gene: str, gene2fam: dict) -> tuple:
    """
    Match a gene to a family using fuzzy string patterns.

    Parameters
    ----------
    gene : str
        Gene identifier to match.
    gene2fam : dict
        Mapping from gene identifiers to family names.

    Returns
    -------
    tuple
        (matched_gene, matched_family) if found, otherwise (None, None).

    Assumptions
    -----------
    Gene identifiers can be matched using regex patterns.
    """
    for key, value in gene2fam.items():
        if fuzzy_match(gene, key):
            return key, value
    return None, None


def get_fam_dic(t: any, gene2fam: dict, new_named_gene2gene_dic: dict) -> dict:
    """
    Map gene families to lists of node names in a tree.

    Parameters
    ----------
    t : any
        Tree object to iterate.
    gene2fam : dict
        Mapping from gene identifiers to family names.
    new_named_gene2gene_dic : dict
        Mapping from node names to gene identifiers.

    Returns
    -------
    dict
        Mapping from family names to node name lists.

    Assumptions
    -----------
    Each tree node name is present in ``new_named_gene2gene_dic``.
    """
    fam_dic = {}
    for i in t:
        gene = new_named_gene2gene_dic[i.name]
        match_gene, match_family = get_matched_value(gene, gene2fam)
        if match_family and match_family in gene2fam.values():
            fam_dic.setdefault(match_family, []).append(i.name)
    return fam_dic


def find_combinations(my_list: list) -> list:
    """
    Generate all pairwise combinations from a list.

    Parameters
    ----------
    my_list : list
        Input list for combination generation.

    Returns
    -------
    list
        List of tuple pairs representing combinations.

    Assumptions
    -----------
    The list length is finite and combinable.
    """
    combinations = []
    for i in range(len(my_list)):
        for j in range(i + 1, len(my_list)):
            combinations.append((my_list[i], my_list[j]))
    return combinations


def get_dup_family_dic(t: any, gene2fam: dict, new_named_gene2gene_dic: dict) -> dict:
    """
    Identify duplication nodes for each gene family.

    Parameters
    ----------
    t : any
        Gene tree object.
    gene2fam : dict
        Mapping from gene identifiers to family labels.
    new_named_gene2gene_dic : dict
        Mapping from node names to gene identifiers.

    Returns
    -------
    dict
        Mapping from family labels to sets of duplication node names.

    Assumptions
    -----------
    Duplication nodes are identified by common ancestors of family members.
    """
    fam_node_dic = {}
    fam_dic = get_fam_dic(t, gene2fam, new_named_gene2gene_dic)
    for k, v in fam_dic.items():
        nodes = set()
        com = find_combinations(v)
        for i in com:
            clade = t.get_common_ancestor(i)
            nodes.add(clade.name)
        fam_node_dic[k] = nodes
    return fam_node_dic


def mapping_sptree(
    t: any,
    sptree: any,
    sp_node_dic: dict,
    gene2fam: dict,
    new_named_gene2gene_dic: dict,
) -> None:
    """
    Map gene family duplication nodes onto a species tree.

    Parameters
    ----------
    t : any
        Gene tree object.
    sptree : any
        Species tree object.
    sp_node_dic : dict
        Mapping from species tree node names to gene family sets.
    gene2fam : dict
        Mapping from gene identifiers to family labels.
    new_named_gene2gene_dic : dict
        Mapping from node names to gene identifiers.

    Returns
    -------
    None

    Assumptions
    -----------
    Species sets can be extracted and matched to species tree nodes.
    """
    fam_node_dic = get_dup_family_dic(t, gene2fam, new_named_gene2gene_dic)
    for k, v in fam_node_dic.items():
        for i in v:
            clade = t & i
            sps = get_species_list(clade)
            uniq_sps = set(sps)
            clade2sptree = sptree.get_common_ancestor(uniq_sps)
            sp_node_dic.setdefault(clade2sptree.name, set()).add(k)


def get_sptree_style(sorted_dict: dict) -> any:
    """
    Create a TreeStyle for species tree rendering with a legend.

    Parameters
    ----------
    sorted_dict : dict
        Mapping from family labels to color and suffix strings.

    Returns
    -------
    any
        Configured TreeStyle object.

    Assumptions
    -----------
    The TreeStyle legend API is available.
    """
    ts = TreeStyle()
    ts.legend_position = 1
    ts.mode = "r"
    ts.scale = 30
    ts.show_border = True
    ts.margin_bottom = 20
    ts.margin_left = 20
    ts.margin_right = 50
    ts.margin_top = 20
    ts.show_leaf_name = True
    ts.extra_branch_line_type = 0
    ts.extra_branch_line_color = "black"
    ts.branch_vertical_margin = -1
    for k, v in sorted_dict.items():
        ts.legend.add_face(
            TextFace(v.split("@")[1] + " " + k, fsize=20, fgcolor=v.split("@")[0]),
            column=0,
        )
    return ts


def create_species_mappings(dict_list: list[dict[str, str]]) -> list[dict[str, str]]:
    """
    Derive species-level mappings from gene-level annotations.

    Parameters
    ----------
    dict_list : list[dict[str, str]]
        List of gene-level mappings: gene2sps, gene2family, gene2order, gene2clade.

    Returns
    -------
    list[dict[str, str]]
        List containing sps2family, sps2order, and sps2clade mappings.

    Assumptions
    -----------
    Gene identifiers are shared across the input dictionaries.
    """
    gene2sps, gene2family, gene2order, gene2clade = dict_list
    sps2family = {}
    sps2order = {}
    sps2clade = {}
    for gene, sps in gene2sps.items():
        if gene in gene2family:
            sps2family[sps] = gene2family[gene]
        if gene in gene2order:
            sps2order[sps] = gene2order[gene]
        if gene in gene2clade:
            sps2clade[sps] = gene2clade[gene]
    return [sps2family, sps2order, sps2clade]


def mark_gene_to_sptree(
    sptree,
    sp_node_dic: dict,
    gene_categories: list,
    sorted_dict: dict,
    gene2sps: dict,
) -> None:
    """
    Mark gene categories on a species tree using colored faces.

    Parameters
    ----------
    sptree : object
        Species tree object to annotate.
    sp_node_dic : dict
        Mapping from species node names to gene family sets.
    gene_categories : list
        List of gene category dictionaries.
    sorted_dict : dict
        Mapping from family labels to color and suffix strings.
    gene2sps : dict
        Mapping from gene identifiers to species names.
    Returns
    -------
    None

    Assumptions
    -----------
    Species names in the tree match those in gene category dictionaries.
    """
    gene_categories_1 = gene_categories.copy()
    gene_categories_1.insert(0, gene2sps)
    gene_categories_2 = create_species_mappings(gene_categories_1)
    color_dicts = generate_color_dict(gene_categories_2)
    faces_added = set()

    def add_face_to_node(node, face, column, position="aligned"):
        if (node, column, position) not in faces_added:
            node.add_face(face, column=column, position=position)
            faces_added.add((color_dict[node.name].split("@")[0], column, position))

    def generate_face_mark(node, species, column, color_dict):
        if (color_dict[node.name].split("@")[0], column, "aligned") not in faces_added:
            if species in color_dict:
                color = color_dict[species].split("@")[ -1]
                face = TextFace(
                    "   ▐" + "  " + color_dict[species].split("@")[0],
                    fgcolor=color,
                    ftype="Arial",
                    fstyle="italic",
                )
                add_face_to_node(node, face, column, position="aligned")
        else:
            color = color_dict[species].split("@")[ -1]
            face = TextFace("   ▐" + "  ", fgcolor=color, ftype="Arial", fstyle="italic")
            add_face_to_node(node, face, column, position="aligned")

    for i in sptree.traverse():
        if not i.is_leaf():
            for k, v in sp_node_dic.items():
                if i.name == k:
                    n = len(v)
                    for index, value in enumerate(sorted(v)):
                        position = (
                            "branch-top" if n == 1 or index < n / 2 else "branch-bottom"
                        )
                        column = index if n == 1 or index < n / 2 else index - n / 2

                        i.add_face(
                            TextFace(
                                sorted_dict[value].split("@")[1],
                                fgcolor=sorted_dict[value].split("@")[0],
                            ),
                            column=column,
                            position=position,
                        )
        else:
            column = 1
            species = i.name
            for color_dict in color_dicts:
                generate_face_mark(i, species, column, color_dict)
                column += 1


# ======================================================
# Section 7: Species Tree Styling
# ======================================================


def rename_sptree(sptree) -> None:
    """
    Apply consistent node styling for species tree visualization.

    Parameters
    ----------
    sptree : object
        Species tree object to be styled.

    Returns
    -------
    None

    Assumptions
    -----------
    The tree uses ETE node style settings.
    """
    for i in sptree.traverse():
        nstyle = NodeStyle()
        nstyle["vt_line_width"] = 1
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_type"] = 0
        nstyle["hz_line_type"] = 0
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        i.set_style(nstyle)


# ======================================================
# Section 8: Main Pipelines (Orchestration)
# ======================================================


def view_main(
    tre_dic,
    sptree,
    gene2new_named_gene_dic,
    voucher2taxa_dic,
    gene_categories,
    tree_style,
    keep_branch,
    new_named_gene2gene_dic,
    gene2fam=None,
    df=None,
    visual: bool = False,
):
    """
    Render annotated gene trees with categorical tip markings.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to file paths.
    sptree : object
        Species tree used for duplication annotation.
    gene2new_named_gene_dic : dict
        Mapping from gene identifiers to renamed identifiers.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.
    gene_categories : list
        List of gene category dictionaries for coloring.
    tree_style : str
        Tree display mode.
    keep_branch : object
        Flag controlling whether to keep original branch lengths.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to gene identifiers.
    gene2fam : dict, optional
        Mapping from gene identifiers to family labels.
    df : object, optional
        DataFrame with heatmap values to render.
    visual : bool, optional
        Whether to display duplication events.

    Returns
    -------
    None

    Assumptions
    -----------
    Input trees are valid and consistent with identifier mappings.
    """
    dir_path = os.path.join(os.getcwd(), "tree_visualizer/")
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path, exist_ok=True)
    color_dicts = generate_color_dict(gene_categories)
    sps_color_dict = get_color_dict(voucher2taxa_dic)
    gene_color_dict = None
    if gene2fam is not None:
        gene_color_dict = get_color_dict(gene2fam)
    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        Phylo_t0 = read_phylo_tree(tre_path)
        Phylo_t0 = rename_input_tre(Phylo_t0, gene2new_named_gene_dic)
        Phylo_t1, ts = get_treestyle(
            Phylo_t0,
            tree_style,
            tre_ID,
            visual,
            species_tree=sptree,
        )
        Phylo_t1.ladderize()
        Phylo_t1.resolve_polytomy(recursive=True)
        Phylo_t1.sort_descendants("support")
        if str(keep_branch) != "1":
            realign_branch_length(Phylo_t1)
            rejust_root_dist(Phylo_t1)

        tips_mark(
            Phylo_t1,
            voucher2taxa_dic,
            color_dicts,
            sps_color_dict,
            tre_ID,
            ts,
            new_named_gene2gene_dic,
            dir_path,
            gene_color_dict,
            df,
        )
        pbar.update(1)
    pbar.close()


def mark_gene_to_sptree_main(
    tre_dic,
    gene_categories,
    sptree,
    gene2fam,
    gene2sps,
    gene2new_named_gene_dic,
    new_named_gene2gene_dic,
    voucher2taxa_dic,
):
    """
    Render a species tree with mapped gene family annotations.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to file paths.
    gene_categories : list
        List of gene category dictionaries.
    sptree : object
        Species tree to annotate.
    gene2fam : dict
        Mapping from gene identifiers to family labels.
    gene2sps : dict
        Mapping from gene identifiers to species names.
    gene2new_named_gene_dic : dict
        Mapping from gene identifiers to renamed identifiers.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to gene identifiers.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.

    Returns
    -------
    None

    Assumptions
    -----------
    Species tree and gene trees share compatible species labels.
    """
    sorted_dict = get_new_sorted_dict(gene2fam)
    num_tre_node(sptree)
    sp_node_dic = {}
    for k, v in tre_dic.items():
        t = read_phylo_tree(v)
        num_tre_node(t)
        t1 = rename_input_tre(t, gene2new_named_gene_dic)
        mapping_sptree(t1, sptree, sp_node_dic, gene2fam, new_named_gene2gene_dic)

    sp2 = rename_input_tre(sptree, voucher2taxa_dic)

    mark_gene_to_sptree(
        sp2,
        sp_node_dic,
        gene_categories,
        sorted_dict,
        gene2sps,
    )

    rename_sptree(sp2)
    ts = get_sptree_style(sorted_dict)
    realign_branch_length(sp2)
    clade_up = sp2.get_children()[0]
    clade_down = sp2.get_children()[1]
    clade_up.dist = 1
    clade_down.dist = get_max_deepth(clade_up)
    sp2.render(file_name="genefamily_map2_sptree.pdf", tree_style=ts)


# ======================================================
# Section 9: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Tree visualization")
    parser.add_argument("--input_GF_list", required=True, help="Gene family list file")
    parser.add_argument("--input_imap", required=True, help="Imap file")
    parser.add_argument("--input_sps_tree", default=None, help="Species tree file (optional)")
    parser.add_argument("--gene_categories", nargs="*", default=[], help="Category mapping files (e.g. genus order)")
    parser.add_argument("--tree_style", default="r", help="Tree style (r=rectangular, c=circular)")
    parser.add_argument("--keep_branch", type=int, default=1, help="Whether to keep branch lengths (1=yes, 0=no)")
    parser.add_argument("--gene2fam", default=None, help="Gene to family mapping file (optional)")
    parser.add_argument("--visual", action="store_true", help="Enable visual output")
    args = parser.parse_args()

    os.makedirs(os.path.join(os.getcwd(), "pdf_result"), exist_ok=True)
    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, _ = (
        gene_id_transfer(args.input_imap)
    )
    gene2fam = read_and_return_dict(args.gene2fam) if args.gene2fam else None
    gene_categories = [read_and_return_dict(f) for f in args.gene_categories]
    tre_dic = read_and_return_dict(args.input_GF_list)
    sptree = PhyloTree(args.input_sps_tree) if args.input_sps_tree else None
    view_main(
        tre_dic,
        sptree,
        gene2new_named_gene_dic,
        voucher2taxa_dic,
        gene_categories,
        args.tree_style,
        args.keep_branch,
        new_named_gene2gene_dic,
        gene2fam=gene2fam,
        visual=args.visual,
    )
