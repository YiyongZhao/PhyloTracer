"""
Monophyletic ortholog filtering and pruning for PhyloTracer gene trees.

This module removes long-branch and insertion artifacts from gene trees,
optionally renders diagnostic PDFs, and writes pruned trees for downstream
orthology analyses in the PhyloTracer pipeline.
"""

import matplotlib.colors as colors
import matplotlib.pyplot as plt
from __init__ import *
from BranchLength_NumericConverter import (
    trans_branch_length,
    write_tree_to_newick,
)
from PyPDF4 import PdfFileReader, PdfFileWriter

# ======================================================
# Section 1: Gene and Tree Renaming Utilities
# ======================================================


def rename_input_single_tre(
    Phylo_t: object,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
) -> object:
    """
    Rename leaf nodes by prefixing species names to gene identifiers.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object whose leaves will be renamed.
    taxa_dic : dict
        Mapping from gene identifiers to species names.
    new_named_gene2gene_dic : dict
        Mapping from renamed node identifiers to original gene names.

    Returns
    -------
    object
        The input tree with updated leaf names.

    Assumptions
    -----------
    Gene identifiers are keys in ``taxa_dic`` and species names are valid
    prefixes for constructing composite labels.
    """
    for node in Phylo_t:
        gene = new_named_gene2gene_dic.get(node.name, node.name)
        species = taxa_dic.get(gene, None)
        if species:
            node.name = species + "_" + node.name
    return Phylo_t


def rename_output_tre(Phylo_t: object, new_named_gene2gene_dic: dict) -> object:
    """
    Restore original gene identifiers after tree processing.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object whose leaves will be renamed back to genes.
    new_named_gene2gene_dic : dict
        Mapping from renamed node identifiers to original gene names.

    Returns
    -------
    object
        The input tree with restored gene identifiers.

    Assumptions
    -----------
    Leaf names follow the ``species_gene`` composite convention.
    """
    for node in Phylo_t:
        parts = node.name.split("_")
        new_str = "_".join(parts[1:]) if len(parts) > 1 else node.name
        node.name = new_named_gene2gene_dic.get(new_str, new_str)
    return Phylo_t


# ======================================================
# Section 2: Copy Number and Visualization Utilities
# ======================================================


def is_multi_copy_tree(Phylo_t: object) -> bool:
    """
    Determine whether a tree contains multiple gene copies per species.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be evaluated.

    Returns
    -------
    bool
        True if the number of leaves exceeds the number of unique species.

    Assumptions
    -----------
    Species identifiers can be parsed by ``get_species_set``.
    """
    leafs = Phylo_t.get_leaf_names()
    uniq_species = get_species_set(Phylo_t)
    if len(leafs) != len(uniq_species):
        return True
    return False


def set_style(
    Phylo_t: object,
    color_dict: dict,
    new_named_gene2gene_dic: dict,
) -> object:
    """
    Apply consistent node styles and colored labels for visualization.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be styled.
    color_dict : dict
        Mapping from gene identifiers to color strings.
    new_named_gene2gene_dic : dict
        Mapping from renamed node identifiers to original gene names.

    Returns
    -------
    object
        The styled tree object.

    Assumptions
    -----------
    Leaf names follow the ``species_gene`` composite convention.
    """
    for node in Phylo_t.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 0
        nstyle["shape"] = "circle"
        nstyle["fgcolor"] = "black"
        node.set_style(nstyle)
        if node.is_leaf():
            parts = node.name.split("_")
            species_name = parts[0]
            new_str = "_".join(parts[1:]) if len(parts) > 1 else node.name
            gene = new_named_gene2gene_dic.get(new_str, new_str)
            color_info = color_dict.get(gene, None)
            if color_info:
                color = color_info.split("*")[-1]
                face = TextFace(
                    species_name + "_" + gene,
                    fgcolor=color,
                    fstyle="italic",
                )
                node.add_face(face, column=0)
    return Phylo_t


def get_color_dict(dictory: dict) -> dict:
    """
    Assign unique colors to taxa for visualization.

    Parameters
    ----------
    dictory : dict
        Mapping whose values define taxa categories to color.

    Returns
    -------
    dict
        Mapping from original keys to composite ``taxon*color`` strings.

    Assumptions
    -----------
    The number of taxa is finite and can be mapped to a colormap.
    """
    colormap = plt.get_cmap("gist_rainbow")
    unique_values = set(dictory.values())
    colors_lst = [
        colors.rgb2hex(colormap(i))
        for i in np.linspace(0, 1, len(unique_values))
    ]
    color_dict = dict(zip(unique_values, colors_lst))
    sps_color_dict = {
        k: v + "*" + color_dict.get(v)
        for k, v in dictory.items()
        if v in color_dict
    }
    return sps_color_dict


# ======================================================
# Section 3: Single-Taxon Clade Utilities
# ======================================================


def get_single_clades(Phylo_t: object, empty_set: set) -> None:
    """
    Collect clades that contain only a single species.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be searched recursively.
    empty_set : set
        Mutable set used to store qualifying clades.

    Returns
    -------
    None

    Assumptions
    -----------
    Species counts can be computed by ``calculate_species_num``.
    """
    if calculate_species_num(Phylo_t) == 1:
        empty_set.add(Phylo_t)
        return
    for child in Phylo_t.get_children():
        get_single_clades(child, empty_set)


def get_node_single_taxa_dict(Phylo_t: object) -> dict:
    """
    Map each species to the list of single-taxon clades containing it.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be searched.

    Returns
    -------
    dict
        Mapping from species names to lists of single-taxon clades.

    Assumptions
    -----------
    Each single-taxon clade contains exactly one species label.
    """

    def get_single_taxa_caldes(empty_set):
        single_taxa_dict = {}
        for clade in empty_set:
            single_taxa = get_species_list(clade)[0]
            single_taxa_dict.setdefault(single_taxa, []).append(clade)
        return single_taxa_dict

    empty_set = set()
    get_single_clades(Phylo_t, empty_set)
    single_taxa_dict = get_single_taxa_caldes(empty_set)
    sorted_dict = dict(
        sorted(
            single_taxa_dict.items(),
            key=lambda item: len(item[1]),
            reverse=True,
        )
    )
    return sorted_dict


def is_single_tree(Phylo_t: object) -> bool:
    """
    Determine whether the tree is composed of single-taxon clades only.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be evaluated.

    Returns
    -------
    bool
        True if each taxon occurs in only one clade, otherwise False.

    Assumptions
    -----------
    Single-taxon clades can be identified by ``get_node_single_taxa_dict``.
    """
    single_taxa_dict = get_node_single_taxa_dict(Phylo_t)

    def check_single(dictionary):
        """
        Check if all values in the dictionary have length one.

        Parameters
        ----------
        dictionary : dict
            Mapping of taxa to lists of clades.

        Returns
        -------
        bool
            True if every taxa maps to a single clade.
        """
        for value in dictionary.values():
            if len(value) != 1:
                return False
        return True

    return check_single(single_taxa_dict)


# ======================================================
# Section 4: PDF Rendering Utilities
# ======================================================


def merge_pdfs_side_by_side(file1: str, file2: str, output_file: str) -> None:
    """
    Merge two single-page PDFs into a single horizontal layout.

    Parameters
    ----------
    file1 : str
        Path to the first PDF file.
    file2 : str
        Path to the second PDF file.
    output_file : str
        Path to the merged output PDF.

    Returns
    -------
    None

    Assumptions
    -----------
    Both input PDFs contain at least one page.
    """
    pdf_writer = PdfFileWriter()
    with open(file1, "rb") as f1, open(file2, "rb") as f2:
        pdf1 = PdfFileReader(f1)
        pdf2 = PdfFileReader(f2)
        page1 = pdf1.getPage(0)
        page2 = pdf2.getPage(0)
        new_page_width = page1.mediaBox.getWidth() + page2.mediaBox.getWidth()
        new_page_height = max(
            page1.mediaBox.getHeight(),
            page2.mediaBox.getHeight(),
        )
        new_page = pdf_writer.addBlankPage(
            new_page_width,
            new_page_height,
        )
        new_page.mergeTranslatedPage(page1, 0, 0)
        new_page.mergeTranslatedPage(page2, page1.mediaBox.getWidth(), 0)
        with open(output_file, "wb") as output:
            pdf_writer.write(output)


# ======================================================
# Section 5: Branch-Length and Insertion Metrics
# ======================================================


def calculate_avg_length(node: object) -> int:
    """
    Compute the average root-to-leaf distance beneath a node.

    Parameters
    ----------
    node : object
        ETE node whose descendant distances will be averaged.

    Returns
    -------
    int
        Average branch length across descendant leaves.

    Assumptions
    -----------
    Distances are defined for all descendant leaves.
    """
    total_length = 0.0
    leaf_count = 0
    for leaf in node.iter_leaves():
        total_length += node.get_distance(leaf)
        leaf_count += 1
    if leaf_count > 0:
        avg_length = total_length / leaf_count
    else:
        avg_length = 0.0
    return avg_length + (node.dist / leaf_count)


def calculate_branch_length_relative_score(node: object) -> int:
    """
    Compute the relative branch-length score against a sister lineage.

    Parameters
    ----------
    node : object
        ETE node for which the relative score is computed.

    Returns
    -------
    int
        Relative branch-length score comparing node and sister lengths.

    Assumptions
    -----------
    The node has a sister lineage or a fallback length of zero is used.
    """
    if not node.is_leaf():
        avg_length = calculate_avg_length(node)
        sister = (
            node.get_sisters()[0]
            if not node.is_root() and node.get_sisters()
            else None
        )
        if sister and not sister.is_leaf():
            sister_length = calculate_avg_length(sister)
        else:
            sister_length = sister.dist if sister else 0.0
        return (avg_length / sister_length) if sister_length != 0 else 0.0
    else:
        branch_length = node.dist
        sister = (
            node.get_sisters()[0]
            if not node.is_root() and node.get_sisters()
            else None
        )
        if sister and not sister.is_leaf():
            sister_length = calculate_avg_length(sister)
        else:
            sister_length = sister.dist if sister else 0.0
        return (branch_length / sister_length) if sister_length != 0 else 0.0


def calculate_insertion_index(node):
    """
    Compute the insertion index of a node relative to its ancestors.

    Parameters
    ----------
    node : object
        ETE node whose insertion index will be computed.

    Returns
    -------
    int
        Minimum insertion ratio across the path to the root.

    Assumptions
    -----------
    Each ancestor has a defined leaf count.
    """
    insertion_index = 1.0
    current_node = node
    while current_node.up:
        current_tips = len(current_node.get_leaves())
        parent_node = current_node.up
        parent_tips = len(parent_node.get_leaves())

        if parent_tips == 0:
            insertion_index = 0.0
            break
        current_insertion_index = len(node.get_leaves()) / current_tips
        insertion_index = min(insertion_index, current_insertion_index)

        current_node = parent_node

    return insertion_index


def get_root_relative_branch_ratio(leaf: object, avg_length: int) -> int:
    """
    Compute the root-relative branch-length ratio for a leaf.

    Parameters
    ----------
    leaf : object
        Leaf node whose branch length is evaluated.
    avg_length : int
        Average branch length across tips.

    Returns
    -------
    int
        Relative ratio of leaf length to the average length.

    Assumptions
    -----------
    Average length is non-zero or handled by caller.
    """
    branch_length = leaf.dist
    return (branch_length - avg_length) / avg_length


def get_sister_relative_branch_ratio(node: object, sister: object) -> int:
    """
    Compute the relative branch-length ratio between a node and its sister.

    Parameters
    ----------
    node : object
        Node whose branch length is evaluated.
    sister : object
        Sister node used for comparison.

    Returns
    -------
    int
        Relative ratio of node length to sister length.

    Assumptions
    -----------
    Missing sister lengths default to zero.
    """
    if node.is_leaf():
        branch_length = node.dist
    else:
        branch_length = calculate_avg_length(node) if node else 0.0

    if sister and not sister.is_leaf():
        sister_length = calculate_avg_length(sister)
    else:
        sister_length = sister.dist if sister else 0.0

    if sister_length != 0:
        return (branch_length - sister_length) / sister_length
    else:
        return 0.0


def calculate_insertion_depth(clade: object, node: object) -> int:
    """
    Compute the topological insertion depth of a node within a clade.

    Parameters
    ----------
    clade : object
        Clade within which the node is evaluated.
    node : object
        Node whose insertion depth is computed.

    Returns
    -------
    int
        Topological distance from clade root to node.

    Assumptions
    -----------
    Topology-only distances are supported by the tree object.
    """
    if clade is None:
        return 0
    return clade.get_distance(node, topology_only=True)


def calculate_insertion_coverage(clade: object, node: object) -> float:
    """
    Compute the insertion coverage ratio of a node within a clade.

    Parameters
    ----------
    clade : object
        Clade containing the node.
    node : object
        Node whose coverage ratio is computed.

    Returns
    -------
    float
        Ratio of node leaf count to clade leaf count.

    Assumptions
    -----------
    Leaf counts are non-zero or a zero ratio is returned.
    """
    if clade is None or node is None:
        return 0.0
    return len(node) / len(clade)


def calculate_insertion(clade: object, node: object) -> int:
    """
    Compute the insertion score as depth divided by coverage.

    Parameters
    ----------
    clade : object
        Clade in which the node is evaluated.
    node : object
        Node to be evaluated.

    Returns
    -------
    int
        Insertion score for the node.

    Assumptions
    -----------
    Coverage is non-zero or a zero score is returned.
    """
    depth = calculate_insertion_depth(clade, node)
    coverage = calculate_insertion_coverage(clade, node)
    return depth / coverage if coverage != 0 else 0


def get_target_clade(clade: object) -> object:
    """
    Identify the nearest ancestor containing exactly two species.

    Parameters
    ----------
    clade : object
        Starting clade for ancestor traversal.

    Returns
    -------
    object
        Ancestor clade containing exactly two species, or None if absent.

    Assumptions
    -----------
    Species counts can be computed by ``get_species_set``.
    """
    node2root = clade.get_ancestors()
    target_clade = None
    for ancestor in node2root:
        if len(get_species_set(ancestor)) > 2:
            break
        elif len(get_species_set(ancestor)) == 2:
            target_clade = ancestor
    return target_clade


def is_ancestor_sister_same(node: object) -> bool:
    """
    Test whether a node's sister matches the sister of its ancestor.

    Parameters
    ----------
    node : object
        Node to be evaluated.

    Returns
    -------
    bool
        True if sister taxa match between node and its ancestor.

    Assumptions
    -----------
    Species identifiers can be extracted by ``get_species_list``.
    """
    sps_list = get_species_list(node)
    sps = sps_list[0] if sps_list else None
    sisters = node.get_sisters()
    sis = sisters[0] if sisters else None
    sis_name_list = get_species_list(sis) if sis else []
    sis_name = sis_name_list[0] if sis_name_list else None
    if sps is None or sis_name is None or sps == sis_name:
        return False
    an = node.up
    if not an or an.is_root():
        return False
    an_sisters = an.get_sisters()
    an_si = an_sisters[0] if an_sisters else None
    an_name_list = get_species_list(an_si) if an_si else []
    an_name = an_name_list[0] if an_name_list else None
    return sis_name == an_name


def get_tips_avg_length(Phylo_t: object) -> int:
    """
    Compute the mean root-to-tip distance across all leaves.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be evaluated.

    Returns
    -------
    int
        Mean root-to-tip distance.

    Assumptions
    -----------
    The tree contains at least one leaf or zero is returned.
    """
    total = sum([Phylo_t.get_distance(leaf) for leaf in Phylo_t])
    count = len(Phylo_t)
    return (total / count) if count else 0


def get_node_avg_length(Phylo_t: object) -> int:
    """
    Compute the mean distance from root to internal nodes.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be evaluated.

    Returns
    -------
    int
        Mean distance from root to internal nodes.

    Assumptions
    -----------
    Internal nodes are identifiable and distances are defined.
    """
    node_length = []
    node_num = 0
    for node in Phylo_t.traverse():
        if not node.is_leaf():
            node_num += 1
            node2root = Phylo_t.get_distance(node)
            node_length.append(node2root)
    return (sum(node_length) / node_num) if node_num else 0


def _resolve_original_gene_name(
    leaf_name: str,
    new_named_gene2gene_dic: dict,
) -> str:
    """
    Resolve original gene identifiers from composite leaf names.

    Parameters
    ----------
    leaf_name : str
        Leaf name possibly prefixed with a species identifier.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.

    Returns
    -------
    str
        Resolved original gene identifier.

    Assumptions
    -----------
    Composite names follow the ``species_gene`` convention.
    """
    parts = leaf_name.split("_")
    candidate = "_".join(parts[1:]) if len(parts) > 1 else leaf_name
    return new_named_gene2gene_dic.get(
        candidate,
        new_named_gene2gene_dic.get(leaf_name, leaf_name),
    )


# ======================================================
# Section 6: Pruning Logic (Core)
# ======================================================


def remove_long_gene(
    Phylo_t: object,
    long_branch_index: int,
    outfile: str,
    tre_ID: str,
    new_named_gene2gene_dic: dict,
) -> object:
    """
    Remove leaves with long branches based on length thresholds.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be pruned.
    long_branch_index : int
        Threshold for detecting long-branch outliers.
    outfile : str
        File handle for recording pruning decisions.
    tre_ID : str
        Identifier of the current tree.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.

    Returns
    -------
    object
        Pruned tree with long-branch leaves removed.

    Assumptions
    -----------
    Root-to-tip and sister-relative ratios capture long-branch behavior.
    """
    Phylo_t1 = Phylo_t.copy()
    remove_gene_set = set()
    tips_avg_length = get_tips_avg_length(Phylo_t1)
    node_avg_length = get_node_avg_length(Phylo_t1)
    while True:
        is_modified = False

        for leaf in Phylo_t1:
            sister = (
                leaf.get_sisters()[0]
                if not leaf.is_root() and leaf.get_sisters()
                else None
            )
            distance = Phylo_t1.get_distance(leaf)
            distance2root_radio = (
                abs(distance / tips_avg_length) if tips_avg_length else 0
            )
            leaf2sister_radio = abs(
                get_sister_relative_branch_ratio(leaf, sister)
            )

            if (
                distance2root_radio > long_branch_index
                or leaf2sister_radio > long_branch_index
            ):
                gene_name = _resolve_original_gene_name(
                    leaf.name,
                    new_named_gene2gene_dic,
                )
                outfile.write(
                    tre_ID
                    + "\t"
                    + "*"
                    + "\t"
                    + gene_name
                    + "\t"
                    + str(distance2root_radio)
                    + "\t"
                    + str(leaf2sister_radio)
                    + "\t"
                    + "\n"
                )
                remove_gene_set.add(leaf.name)
                is_modified = True
            else:
                gene_name = _resolve_original_gene_name(
                    leaf.name,
                    new_named_gene2gene_dic,
                )
                outfile.write(
                    tre_ID
                    + "\t"
                    + "\t"
                    + "\t"
                    + gene_name
                    + "\t"
                    + str(distance2root_radio)
                    + "\t"
                    + str(leaf2sister_radio)
                    + "\t"
                    + "\n"
                )

        if not is_modified:
            break

        total_leafs_set = set(Phylo_t1.get_leaf_names())
        diff = total_leafs_set - remove_gene_set
        Phylo_t1.prune(diff, preserve_branch_length=True)

    return Phylo_t1


def remove_insert_gene(
    Phylo_t: object,
    long_branch_index: int,
    outfile: str,
    tre_ID: str,
    new_named_gene2gene_dic: dict,
) -> object:
    """
    Remove inserted genes based on insertion indices and topology.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be pruned.
    long_branch_index : int
        Threshold for insertion-related pruning.
    outfile : str
        File handle for recording pruning decisions.
    tre_ID : str
        Identifier of the current tree.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.

    Returns
    -------
    object
        Pruned tree with insertion artifacts removed.

    Assumptions
    -----------
    Insertion depth and coverage metrics capture abnormal placements.
    """
    Phylo_t1 = Phylo_t.copy()
    taxa_clade = set()
    get_single_clades(Phylo_t1, taxa_clade)
    remove_gene_set = set()

    for clade in taxa_clade:
        temp_name = "_".join(clade.name.split("_")[1:])
        if is_ancestor_sister_same(clade):
            target_clade = get_target_clade(clade)
            if target_clade is None:
                continue
            else:
                index_num = calculate_insertion_depth(target_clade, clade)
                index_over = calculate_insertion_coverage(target_clade, clade)
                insert = calculate_insertion(target_clade, clade)
                if clade.is_leaf():
                    gene_name = _resolve_original_gene_name(
                        clade.name,
                        new_named_gene2gene_dic,
                    )
                    outfile.write(
                        tre_ID
                        + "\t"
                        + "@"
                        + "\t"
                        + gene_name
                        + "\t"
                        + str(index_num)
                        + "\t"
                        + str(index_over)
                        + "\t"
                        + str(insert)
                        + "\n"
                    )
                else:
                    outfile.write(
                        tre_ID
                        + "\t"
                        + "@"
                        + "\t"
                        + clade.name
                        + "\t"
                        + str(index_num)
                        + "\t"
                        + str(index_over)
                        + "\t"
                        + str(insert)
                        + "\n"
                    )
                remove_gene_set.add(clade.name)
    total_leafs_set = set(Phylo_t1.get_leaf_names())
    diff = total_leafs_set - remove_gene_set
    Phylo_t1.prune(diff, preserve_branch_length=True)
    return Phylo_t1


def prune_single(Phylo_t):
    """
    Prune single-taxon clades to reduce redundant taxa copies.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be pruned.

    Returns
    -------
    object
        Pruned tree with single-taxon clades removed.

    Assumptions
    -----------
    Single-taxon clades are identifiable by ``get_node_single_taxa_dict``.
    """
    rm_list = []
    single_taxa_dict = get_node_single_taxa_dict(Phylo_t)
    for k, v in single_taxa_dict.items():
        if len(v) == 1:
            pass
        elif len(v) == 2:
            if len(v[0]) > len(v[1]):
                leafs = v[1].get_leaf_names()
                total_leafs = Phylo_t.get_leaf_names()
                diff = [a for a in total_leafs if a not in set(leafs)]
                Phylo_t.prune(diff, preserve_branch_length=True)

            else:
                leafs = v[0].get_leaf_names()
                total_leafs = Phylo_t.get_leaf_names()
                diff = [a for a in total_leafs if a not in set(leafs)]
                Phylo_t.prune(diff, preserve_branch_length=True)
        else:
            result = []
            for i in v:
                insertion_index = calculate_insertion_index(i)
                result.append(insertion_index)
            min_score = min(result)
            d = []
            for j in v:
                insertion_index1 = calculate_insertion_index(j)
                if insertion_index1 == min_score:
                    branch_length_relative_score = (
                        calculate_branch_length_relative_score(j)
                    )
                    d.append(branch_length_relative_score)
            min_score1 = min(d)
            for h in v:
                branch_length_relative_score1 = (
                    calculate_branch_length_relative_score(h)
                )
                if branch_length_relative_score1 == min_score1:
                    rm_list.append(h.name)

    for i in rm_list:
        clade = Phylo_t & i
        if clade.is_leaf():
            clade.delete()
        else:
            leafs = clade.get_leaf_names()
            total_leafs = Phylo_t.get_leaf_names()
            diff = [a for a in total_leafs if a not in set(leafs)]
            Phylo_t.prune(diff, preserve_branch_length=True)


# ======================================================
# Section 7: Visualization Outputs
# ======================================================


def generate_pdf_before(tre_ID, Phylo_t) -> None:
    """
    Render a PDF snapshot of the tree before pruning.

    Parameters
    ----------
    tre_ID : str
        Tree identifier used in the output filename.
    Phylo_t : object
        ETE tree object to render.

    Returns
    -------
    None

    Assumptions
    -----------
    The ETE rendering backend is available.
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(tre_ID + "_before", fsize=10), column=0)
    Phylo_t.render(file_name=tre_ID + "_before.pdf", tree_style=ts)


def generate_pdf_after(tre_ID, Phylo_t) -> None:
    """
    Render a PDF snapshot of the tree after pruning.

    Parameters
    ----------
    tre_ID : str
        Tree identifier used in the output filename.
    Phylo_t : object
        ETE tree object to render.

    Returns
    -------
    None

    Assumptions
    -----------
    The ETE rendering backend is available.
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(tre_ID + "_after", fsize=10), column=0)
    Phylo_t.render(file_name=tre_ID + "_after.pdf", tree_style=ts)


# ======================================================
# Section 8: Main Pipeline (Orchestration)
# ======================================================


def prune_main_Mono(
    tre_dic: dict,
    taxa_dic: dict,
    long_branch_index: int,
    insert_branch_index: int,
    new_named_gene2gene_dic: dict,
    gene2new_named_gene_dic: dict,
    visual: bool = False,
) -> None:
    """
    Run mono-copy ortholog pruning with optional visualization.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to tree file paths.
    taxa_dic : dict
        Mapping from gene names to species names.
    long_branch_index : int
        Threshold for identifying long branches.
    insert_branch_index : int
        Threshold for identifying insertion artifacts.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.
    gene2new_named_gene_dic : dict
        Mapping from original gene names to renamed identifiers.
    visual : bool, optional
        Whether to generate and merge diagnostic PDFs.

    Returns
    -------
    None

    Assumptions
    -----------
    Input trees are valid Newick files and identifiers are mappable.
    """
    color_dic = get_color_dict(taxa_dic)
    dir_path1 = os.path.join(os.getcwd(), "orthofilter_mono/pruned_tree/")
    if os.path.exists(dir_path1):
        shutil.rmtree(dir_path1)
    os.makedirs(dir_path1)
    if visual:
        dir_path2 = os.path.join(os.getcwd(), "orthofilter_mono/pruned_tree_pdf/")
        if os.path.exists(dir_path2):
            shutil.rmtree(dir_path2)
        os.makedirs(dir_path2)
    dir_path3 = os.path.join(os.getcwd(), "orthofilter_mono/long_branch_gene/")
    if os.path.exists(dir_path3):
        shutil.rmtree(dir_path3)
    os.makedirs(dir_path3)

    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t0 = Tree(tre_path)
        t0.ladderize()
        t0.resolve_polytomy(recursive=True)
        t0.sort_descendants("support")
        t = rename_input_tre(t0, gene2new_named_gene_dic)
        num_tre_node(t)

        if is_multi_copy_tree(t):
            o = open(
                os.path.join(dir_path3, tre_ID + "_delete_gene.txt"),
                "w",
            )
            o.write(
                "tre_ID"
                + "\t"
                + "long_branch_label"
                + "\t"
                + "gene"
                + "\t"
                + "root_relative_branch_ratio"
                + "\t"
                + "sister_relative_branch_ratio"
                + "\n"
            )
            rename_input_single_tre(t, taxa_dic, new_named_gene2gene_dic)

            if visual:
                set_style(t, color_dic, new_named_gene2gene_dic)
                generate_pdf_before(tre_ID, t)
            t1 = t
            # t1 = remove_long_gene(t, long_branch_index, o, tre_ID,new_named_gene2gene_dic)
            o.write("\n")
            if len(get_species_set(t1)) != 1:
                o.write(
                    "tre_ID"
                    + "\t"
                    + "insert_branch_label"
                    + "\t"
                    + "gene"
                    + "\t"
                    + "insertion_depth"
                    + "\t"
                    + "insertion_coverage"
                    + "\t"
                    + "calculate_insertion"
                    + "\n"
                )
                t2 = remove_insert_gene(
                    t1,
                    insert_branch_index,
                    o,
                    tre_ID,
                    new_named_gene2gene_dic,
                )
                o.close()
                if visual:
                    generate_pdf_after(tre_ID, t2)
                    merge_pdfs_side_by_side(
                        tre_ID + "_before.pdf",
                        tre_ID + "_after.pdf",
                        os.path.join(dir_path2, tre_ID + ".pdf"),
                    )
                    os.remove(tre_ID + "_before.pdf")
                    os.remove(tre_ID + "_after.pdf")
                t3 = rename_output_tre(t2, new_named_gene2gene_dic)
                tree_str = trans_branch_length(t3)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)
            else:
                if visual:
                    generate_pdf_after(tre_ID, t1)
                    merge_pdfs_side_by_side(
                        tre_ID + "_before.pdf",
                        tre_ID + "_after.pdf",
                        os.path.join(dir_path2, tre_ID + ".pdf"),
                    )
                    os.remove(tre_ID + "_before.pdf")
                    os.remove(tre_ID + "_after.pdf")
                t3 = rename_output_tre(t1, new_named_gene2gene_dic)
                tree_str = trans_branch_length(t3)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)

        else:
            o = open(
                os.path.join(dir_path3, tre_ID + "_delete_gene.txt"),
                "w",
            )
            o.write(
                "tre_ID"
                + "\t"
                + "long_branch_label"
                + "\t"
                + "gene"
                + "\t"
                + "root_relative_branch_ratio"
                + "\t"
                + "sister_relative_branch_ratio"
                + "\n"
            )
            rename_input_single_tre(t, taxa_dic, new_named_gene2gene_dic)
            if visual:
                set_style(t, color_dic, new_named_gene2gene_dic)
                generate_pdf_before(tre_ID, t)
            t1 = t
            # t1 = remove_long_gene(t, long_branch_index, o, tre_ID,new_named_gene2gene_dic)
            o.write("\n")
            if len(get_species_set(t1)) != 1:
                o.write(
                    "tre_ID"
                    + "\t"
                    + "insert_branch_label"
                    + "\t"
                    + "gene"
                    + "\t"
                    + "insertion_depth"
                    + "\t"
                    + "insertion_coverage"
                    + "\t"
                    + "calculate_insertion"
                    + "\n"
                )
                t2 = remove_insert_gene(
                    t1,
                    insert_branch_index,
                    o,
                    tre_ID,
                    new_named_gene2gene_dic,
                )

                while not is_single_tree(t2):
                    prune_single(t2)
                    pass
                if visual:
                    generate_pdf_after(tre_ID, t2)
                    merge_pdfs_side_by_side(
                        tre_ID + "_before.pdf",
                        tre_ID + "_after.pdf",
                        os.path.join(dir_path2, tre_ID + ".pdf"),
                    )
                    os.remove(tre_ID + "_before.pdf")
                    os.remove(tre_ID + "_after.pdf")

                t3 = rename_output_tre(t2, new_named_gene2gene_dic)
                tree_str = trans_branch_length(t3)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)
            else:
                if visual:
                    generate_pdf_after(tre_ID, t1)
                    merge_pdfs_side_by_side(
                        tre_ID + "_before.pdf",
                        tre_ID + "_after.pdf",
                        os.path.join(dir_path2, tre_ID + ".pdf"),
                    )
                    os.remove(tre_ID + "_before.pdf")
                    os.remove(tre_ID + "_after.pdf")
                t3 = rename_output_tre(t1, new_named_gene2gene_dic)
                tree_str = trans_branch_length(t3)
                write_tree_to_newick(tree_str, tre_ID, dir_path1)

            o.close()
        pbar.update(1)
    pbar.close()


# ======================================================
# Section 9: Summary Wrapper and CLI Entry Point
# ======================================================


def prune_mono_copy_trees(
    tre_dic: dict,
    taxa_dic: dict,
    long_branch_index: int,
    insert_branch_index: int,
    new_named_gene2gene_dic: dict,
    gene2new_named_gene_dic: dict,
    visual: bool = False,
) -> None:
    """
    Wrapper for mono-copy pruning to preserve legacy entry points.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to tree file paths.
    taxa_dic : dict
        Mapping from gene names to species names.
    long_branch_index : int
        Threshold for identifying long branches.
    insert_branch_index : int
        Threshold for identifying insertion artifacts.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.
    gene2new_named_gene_dic : dict
        Mapping from original gene names to renamed identifiers.
    visual : bool, optional
        Whether to generate diagnostic PDFs.

    Returns
    -------
    None

    Assumptions
    -----------
    This wrapper preserves the existing function signature and behavior.
    """
    return prune_main_Mono(
        tre_dic,
        taxa_dic,
        long_branch_index,
        insert_branch_index,
        new_named_gene2gene_dic,
        gene2new_named_gene_dic,
        visual,
    )


if __name__ == "__main__":
    os.makedirs(os.path.join(os.getcwd(), "pruned_tree"))
    taxa_dic = read_and_return_dict("taxa.txt")
    tre_dic = read_and_return_dict("100_nosingle_GF_list.txt")
    long_brancch_index = 5
    insert_branch_index = 5
    prune_main_Mono(tre_dic, taxa_dic, long_brancch_index, insert_branch_index)
