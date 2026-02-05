"""
Monophyletic ortholog pruning based on insertion topology for PhyloTracer.

This module removes insertion artifacts from gene trees using topology-based
insertion depth and coverage metrics. It optionally renders diagnostic PDFs
and writes pruned trees for downstream orthology inference.

Long-branch filtering is intentionally excluded and handled in a separate module.
"""

# ======================================================
# Imports
# ======================================================

from __init__ import *

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from PyPDF4 import PdfFileReader, PdfFileWriter
from BranchLength_NumericConverter import (
    trans_branch_length,
    write_tree_to_newick,
)

# ======================================================
# Section 1: Naming & Tree Utilities
# ======================================================


def rename_input_single_tre(
    Phylo_t: object,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
) -> object:
    """
    Prefix species names to gene identifiers on leaf nodes.

    Parameters
    ----------
    Phylo_t : object
        ETE gene tree whose leaves will be renamed.
    taxa_dic : dict
        Mapping from gene identifiers to species names.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.

    Returns
    -------
    object
        Tree with leaf names updated as ``species_gene``.

    Assumptions
    -----------
    Gene identifiers can be mapped to species using ``taxa_dic``.
    """
    for node in Phylo_t:
        gene = new_named_gene2gene_dic.get(node.name, node.name)
        species = taxa_dic.get(gene)
        if species:
            node.name = f"{species}_{node.name}"
    return Phylo_t


def rename_output_tre(
    Phylo_t: object,
    new_named_gene2gene_dic: dict,
) -> object:
    """
    Restore original gene identifiers from composite leaf names.

    Parameters
    ----------
    Phylo_t : object
        ETE gene tree whose leaves will be renamed.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.

    Returns
    -------
    object
        Tree with restored gene identifiers.

    Assumptions
    -----------
    Leaf names follow the ``species_gene`` composite convention.
    """
    for node in Phylo_t:
        parts = node.name.split("_")
        gene = "_".join(parts[1:]) if len(parts) > 1 else node.name
        node.name = new_named_gene2gene_dic.get(gene, gene)
    return Phylo_t


def _resolve_original_gene_name(
    leaf_name: str,
    new_named_gene2gene_dic: dict,
) -> str:
    """
    Resolve original gene identifier from a composite leaf name.

    Parameters
    ----------
    leaf_name : str
        Leaf name encoded as ``species_gene`` or a raw gene ID.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.

    Returns
    -------
    str
        Original gene identifier.

    Assumptions
    -----------
    Composite names use underscores to separate species and gene parts.
    """
    parts = leaf_name.split("_")
    candidate = "_".join(parts[1:]) if len(parts) > 1 else leaf_name
    return new_named_gene2gene_dic.get(candidate, leaf_name)


# ======================================================
# Section 2: Copy Number & Taxon Structure
# ======================================================


def is_multi_copy_tree(Phylo_t: object) -> bool:
    """
    Determine whether a tree contains multiple copies per species.

    Parameters
    ----------
    Phylo_t : object
        ETE gene tree to evaluate.

    Returns
    -------
    bool
        True if leaf count differs from the number of unique species.

    Assumptions
    -----------
    Species identifiers can be parsed by ``get_species_set``.
    """
    return len(Phylo_t.get_leaf_names()) != len(get_species_set(Phylo_t))


def get_single_clades(Phylo_t: object, result_set: set) -> None:
    """
    Recursively collect clades containing a single species.

    Parameters
    ----------
    Phylo_t : object
        ETE tree or subtree to traverse.
    result_set : set
        Mutable set to collect single-taxon clades.

    Returns
    -------
    None

    Assumptions
    -----------
    Species counts are computed by ``get_species_set``.
    """
    if len(get_species_set(Phylo_t)) == 1:
        result_set.add(Phylo_t)
        return
    for child in Phylo_t.get_children():
        get_single_clades(child, result_set)


def get_node_single_taxa_dict(Phylo_t: object) -> dict:
    """
    Map each species to its corresponding single-taxon clades.

    Parameters
    ----------
    Phylo_t : object
        ETE gene tree to analyze.

    Returns
    -------
    dict
        Mapping from species to lists of single-taxon clades.

    Assumptions
    -----------
    Each single-taxon clade contains exactly one species label.
    """
    result = {}
    clades = set()
    get_single_clades(Phylo_t, clades)

    for clade in clades:
        species = get_species_list(clade)[0]
        result.setdefault(species, []).append(clade)

    return dict(sorted(result.items(), key=lambda x: len(x[1]), reverse=True))


def is_single_tree(Phylo_t: object) -> bool:
    """
    Check whether each species appears in only one clade.

    Parameters
    ----------
    Phylo_t : object
        ETE gene tree to evaluate.

    Returns
    -------
    bool
        True if each species maps to exactly one clade.

    Assumptions
    -----------
    Single-taxon clades are identified by ``get_node_single_taxa_dict``.
    """
    return all(len(v) == 1 for v in get_node_single_taxa_dict(Phylo_t).values())


# ======================================================
# Section 3: Insertion Metrics (Topology-Based)
# ======================================================


def calculate_insertion_index(node: object) -> float:
    """
    Compute minimum insertion ratio along the path to root.

    Parameters
    ----------
    node : object
        ETE node for which the insertion ratio is computed.

    Returns
    -------
    float
        Minimum ratio of node leaf count to ancestor leaf count.

    Assumptions
    -----------
    Node leaf counts are non-zero along the ancestor path.
    """
    index = 1.0
    current = node

    while current.up:
        parent = current.up
        if len(parent) == 0:
            return 0.0
        ratio = len(node) / len(current)
        index = min(index, ratio)
        current = parent

    return index


def calculate_insertion_depth(clade: object, node: object) -> int:
    """
    Compute topological distance between clade root and node.

    Parameters
    ----------
    clade : object
        Clade used as the reference root.
    node : object
        Node whose depth is measured within the clade.

    Returns
    -------
    int
        Topological distance from clade to node.

    Assumptions
    -----------
    Topology-only distances are supported by the tree object.
    """
    return clade.get_distance(node, topology_only=True) if clade else 0


def calculate_insertion_coverage(clade: object, node: object) -> float:
    """
    Compute coverage ratio of node tips relative to clade tips.

    Parameters
    ----------
    clade : object
        Clade used as the reference set of tips.
    node : object
        Node whose tip coverage is measured.

    Returns
    -------
    float
        Ratio of node tips to clade tips.

    Assumptions
    -----------
    Clade size is non-zero when defined.
    """
    return len(node) / len(clade) if clade and len(clade) else 0.0


def calculate_insertion(clade: object, node: object) -> float:
    """
    Compute insertion score as depth divided by coverage.

    Parameters
    ----------
    clade : object
        Reference clade for depth calculation.
    node : object
        Node to score.

    Returns
    -------
    float
        Insertion score for the node.

    Assumptions
    -----------
    Coverage values of zero yield a score of 0.0.
    """
    coverage = calculate_insertion_coverage(clade, node)
    return calculate_insertion_depth(clade, node) / coverage if coverage != 0 else 0.0


def get_target_clade(clade: object) -> object:
    """
    Identify nearest ancestor containing exactly two species.

    Parameters
    ----------
    clade : object
        Starting clade for ancestor traversal.

    Returns
    -------
    object
        Target ancestor clade or None if not found.

    Assumptions
    -----------
    Species counts are computed by ``get_species_set``.
    """
    for ancestor in clade.get_ancestors():
        sp_num = len(get_species_set(ancestor))
        if sp_num > 2:
            break
        if sp_num == 2:
            return ancestor
    return None


def is_ancestor_sister_same(node: object) -> bool:
    """
    Test whether node and its ancestor share the same sister taxon.

    Parameters
    ----------
    node : object
        Node to evaluate.

    Returns
    -------
    bool
        True if the sister taxa match across the node and its ancestor.

    Assumptions
    -----------
    Species identifiers can be extracted by ``get_species_list``.
    """
    if not node.up or node.up.is_root():
        return False

    sps = get_species_list(node)
    if not sps:
        return False

    sister = node.get_sisters()[0] if node.get_sisters() else None
    ancestor_sister = node.up.get_sisters()[0] if node.up.get_sisters() else None

    if not sister or not ancestor_sister:
        return False

    return get_species_list(sister)[0] == get_species_list(ancestor_sister)[0]


# ======================================================
# Section 4: Pruning Core (Insertion-based)
# ======================================================


def remove_insert_gene(
    Phylo_t: object,
    inserted_depth: int,
    outfile,
    tre_ID: str,
    new_named_gene2gene_dic: dict,
) -> object:
    """
    Remove insertion artifacts based on insertion metrics.

    Parameters
    ----------
    Phylo_t : object
        ETE gene tree to be pruned.
    inserted_depth : int
        Depth threshold for removal.
    outfile : object
        File handle for logging removed nodes.
    tre_ID : str
        Tree identifier used in logs.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.

    Returns
    -------
    object
        Pruned gene tree.

    Assumptions
    -----------
    Insertion depth and coverage capture abnormal placements.
    """
    t = Phylo_t.copy()
    single_clades = set()
    get_single_clades(t, single_clades)

    remove_set = set()

    for clade in single_clades:
        if not is_ancestor_sister_same(clade):
            continue

        target = get_target_clade(clade)
        if not target:
            continue

        depth = calculate_insertion_depth(target, clade)
        coverage = calculate_insertion_coverage(target, clade)
        score = calculate_insertion(target, clade)

        gene = _resolve_original_gene_name(clade.name, new_named_gene2gene_dic)

        outfile.write(f"{tre_ID}\t@\t{gene}\t{depth}\t{coverage}\t{score}\n")

        if depth >= inserted_depth:
            remove_set.add(clade.name)

    keep = set(t.get_leaf_names()) - remove_set
    t.prune(keep, preserve_branch_length=True)
    return t


def prune_single(Phylo_t: object) -> None:
    """
    Resolve remaining redundant single-taxon clades.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object to be pruned in place.

    Returns
    -------
    None

    Assumptions
    -----------
    Insertion indices prioritize clades for removal.
    """
    rm_list = []
    single_taxa = get_node_single_taxa_dict(Phylo_t)

    for _, clades in single_taxa.items():
        if len(clades) <= 1:
            continue

        scores = [calculate_insertion_index(c) for c in clades]
        min_score = min(scores)

        for c in clades:
            if calculate_insertion_index(c) == min_score:
                rm_list.append(c.name)

    for name in rm_list:
        clade = Phylo_t & name
        if clade.is_leaf():
            clade.delete()
        else:
            keep = set(Phylo_t.get_leaf_names()) - set(clade.get_leaf_names())
            Phylo_t.prune(keep, preserve_branch_length=True)


# ======================================================
# Section 5: Visualization Utilities (Optional)
# ======================================================


def get_color_dict(taxa_dic: dict) -> dict:
    """
    Assign colors to taxa for visualization.

    Parameters
    ----------
    taxa_dic : dict
        Mapping from gene identifiers to taxa labels.

    Returns
    -------
    dict
        Mapping from identifiers to ``taxon*color`` strings.

    Assumptions
    -----------
    A rainbow colormap is used to assign unique taxa colors.
    """
    cmap = plt.get_cmap("gist_rainbow")
    uniq = list(set(taxa_dic.values()))
    color_list = [colors.rgb2hex(cmap(i)) for i in np.linspace(0, 1, len(uniq))]
    return {k: f"{v}*{dict(zip(uniq, color_list))[v]}" for k, v in taxa_dic.items()}


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
                face = TextFace(species_name + "_" + gene, fgcolor=color, fstyle="italic")
                node.add_face(face, column=0)
    return Phylo_t


def generate_pdf(tre_ID: str, Phylo_t: object, tag: str) -> None:
    """
    Render a tree snapshot as a PDF file.

    Parameters
    ----------
    tre_ID : str
        Tree identifier for output naming.
    Phylo_t : object
        Tree to render.
    tag : str
        Tag appended to the output filename.

    Returns
    -------
    None

    Assumptions
    -----------
    ETE rendering is available in the runtime environment.
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(f"{tre_ID}_{tag}", fsize=10), 0)
    Phylo_t.render(f"{tre_ID}_{tag}.pdf", tree_style=ts)


def merge_pdfs_side_by_side(f1: str, f2: str, out: str) -> None:
    """
    Merge two single-page PDFs horizontally.

    Parameters
    ----------
    f1 : str
        Path to the first PDF.
    f2 : str
        Path to the second PDF.
    out : str
        Path to the output merged PDF.

    Returns
    -------
    None

    Assumptions
    -----------
    Input PDFs contain at least one page.
    """
    writer = PdfFileWriter()
    with open(f1, "rb") as a, open(f2, "rb") as b:
        p1 = PdfFileReader(a).getPage(0)
        p2 = PdfFileReader(b).getPage(0)
        page = writer.addBlankPage(
            p1.mediaBox.getWidth() + p2.mediaBox.getWidth(),
            max(p1.mediaBox.getHeight(), p2.mediaBox.getHeight()),
        )
        page.mergeTranslatedPage(p1, 0, 0)
        page.mergeTranslatedPage(p2, p1.mediaBox.getWidth(), 0)
        with open(out, "wb") as o:
            writer.write(o)


# ======================================================
# Section 6: Main Pipeline (Orchestration)
# ======================================================


def prune_main_Mono(
    tre_dic: dict,
    taxa_dic: dict,
    inserted_depth: int,
    new_named_gene2gene_dic: dict,
    gene2new_named_gene_dic: dict,
    visual: bool = False,
) -> None:
    """
    Run the insertion-based mono-copy pruning pipeline.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to file paths.
    taxa_dic : dict
        Mapping from gene identifiers to taxa labels.
    inserted_depth : int
        Depth threshold for insertion pruning.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene names.
    gene2new_named_gene_dic : dict
        Mapping from original gene names to renamed identifiers.
    visual : bool, optional
        Whether to generate before/after PDF visualizations.

    Returns
    -------
    None

    Assumptions
    -----------
    Input trees are valid and compatible with identifier mappings.
    """
    color_dic = get_color_dict(taxa_dic)

    out_tree_dir = "orthofilter_mono/pruned_tree"
    out_log_dir = "orthofilter_mono/insert_gene"
    out_visual_dir = "orthofilter_mono/visual"
    shutil.rmtree(out_tree_dir, ignore_errors=True)
    shutil.rmtree(out_log_dir, ignore_errors=True)
    shutil.rmtree(out_visual_dir, ignore_errors=True)
    os.makedirs(out_tree_dir)
    os.makedirs(out_log_dir)
    os.makedirs(out_visual_dir)

    pbar = tqdm(total=len(tre_dic), desc="Processing trees", unit="tree")
    for tre_ID, tre_path in tre_dic.items():
        pbar.set_description(f"Processing {tre_ID}")
        t0 = Tree(tre_path)
        t0.ladderize()
        t0.resolve_polytomy(recursive=True)
        t0.sort_descendants("support")


        t = rename_input_tre(t0, gene2new_named_gene_dic)
        num_tre_node(t)
        rename_input_single_tre(t, taxa_dic, new_named_gene2gene_dic)
        log = open(os.path.join(out_log_dir, f"{tre_ID}_insert_gene.txt"), "w")
        log.write("tre_ID\tlabel\tgene\tdepth\tcoverage\tinsertion\n")

        if visual:
            set_style(t, color_dic, new_named_gene2gene_dic)
            generate_pdf(tre_ID, t, "before")

        if len(get_species_set(t)) > 1:
            t1 = remove_insert_gene(
                t,
                inserted_depth,
                log,
                tre_ID,
                new_named_gene2gene_dic,
            )

        if visual:
            generate_pdf(tre_ID, t1, "after")
            merge_pdfs_side_by_side(
                f"{tre_ID}_before.pdf",
                f"{tre_ID}_after.pdf",
                os.path.join(out_visual_dir, f"{tre_ID}.pdf"),
            )
            os.remove(f"{tre_ID}_before.pdf")
            os.remove(f"{tre_ID}_after.pdf")

        t2 = rename_output_tre(t1, new_named_gene2gene_dic)
        tree_str = trans_branch_length(t2)
        write_tree_to_newick(tree_str, tre_ID, out_tree_dir)

        log.close()
        pbar.update(1)


# ======================================================
# Section 7: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    taxa_dic = read_and_return_dict("taxa.txt")
    tre_dic = read_and_return_dict("100_nosingle_GF_list.txt")
    prune_main_Mono(
        tre_dic,
        taxa_dic,
        inserted_depth=5,
        new_named_gene2gene_dic={},
        gene2new_named_gene_dic={},
    )
