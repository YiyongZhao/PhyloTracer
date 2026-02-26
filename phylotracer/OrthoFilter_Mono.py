"""
Monophyletic ortholog pruning guided by dominant-lineage purity in PhyloTracer.

This module detects dominant Brassicaceae lineages in gene trees and removes
alien tips using species-tree distance, coverage, and insertion depth metrics.
It optionally renders diagnostic PDFs and writes pruned trees for downstream
orthology inference.

Long-branch filtering is intentionally excluded and handled in a separate module.
"""

import os
import shutil
from collections import Counter

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from ete3 import PhyloTree, Tree
try:
    from ete3 import NodeStyle, TextFace, TreeStyle
except ImportError:
    NodeStyle = None
    TextFace = None
    TreeStyle = None
try:
    from pypdf import PdfReader, PdfWriter
except ImportError:
    PdfReader = None
    PdfWriter = None
from tqdm import tqdm

from phylotracer import (
    num_tre_node,
    read_and_return_dict,
    rename_input_tre,
)
from phylotracer.BranchLength_NumericConverter import (
    trans_branch_length,
    write_tree_to_newick,
)

# --------------------------
# 1. Naming and Tree Utilities
# --------------------------


def rename_input_single_tre(
    Phylo_t: object,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
) -> object:
    """Prefix species names to gene identifiers on leaf nodes.

    Args:
        Phylo_t (object): ETE gene tree whose leaves will be renamed.
        taxa_dic (dict): Mapping from gene identifiers to species names.
        new_named_gene2gene_dic (dict): Mapping from renamed identifiers to original gene names.

    Returns:
        object: Tree with leaf names updated as ``species_gene``.

    Assumptions:
        Gene identifiers can be mapped to species using ``taxa_dic``.
    """
    for node in Phylo_t:
        gene = new_named_gene2gene_dic.get(node.name, node.name)
        lineage = taxa_dic.get(gene)
        if lineage:
            node.name = f"{lineage}_{node.name}"
    return Phylo_t


def rename_output_tre(
    Phylo_t: object,
    new_named_gene2gene_dic: dict,
) -> object:
    """Restore original gene identifiers from composite leaf names.

    Args:
        Phylo_t (object): ETE gene tree whose leaves will be renamed.
        new_named_gene2gene_dic (dict): Mapping from renamed identifiers to original gene names.

    Returns:
        object: Tree with restored gene identifiers.

    Assumptions:
        Leaf names follow the ``species|gene`` composite convention.
    """
    for node in Phylo_t:
        parts = node.name.split("_")
        gene = "_".join(parts[1:])
        node.name = new_named_gene2gene_dic.get(gene, gene)
    return Phylo_t


# --------------------------
# 2. Dominant-Lineage Metrics
# --------------------------


def parse_leaf_components(leaf_name: str) -> tuple:
    """Parse clade and renamed-gene components from a leaf name.

    Args:
        leaf_name (str): Leaf name encoded as ``clade_renamed_gene``.

    Returns:
        tuple: (clade_label, renamed_gene_id). If no clade prefix is detected,
            returns ("", original_name).

    Assumptions:
        The first underscore separates clade prefix from renamed gene ID.
    """
    parts = leaf_name.split("_")
    if len(parts) < 2:
        # Fallback path when no clade prefix exists.
        return "", leaf_name

    gene_id = "_".join(parts[1:])
    clade = parts[0]
    return clade, gene_id


def get_leaf_clade(
    leaf_name: str,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
) -> str:
    """Resolve the clade label for a leaf node.

    Args:
        leaf_name (str): Leaf name encoded as ``clade_renamed_gene``.
        taxa_dic (dict): Mapping from original gene IDs to clade labels.
        new_named_gene2gene_dic (dict): Mapping from renamed identifiers to original gene names.

    Returns:
        str: Clade label for the leaf.

    Assumptions:
        Leaf names were prefixed with clade labels using ``rename_input_single_tre``.
    """
    clade, gene_id = parse_leaf_components(leaf_name)
    if clade:
        return clade

    # Fallback path: map renamed ID back to original gene ID.
    orig = new_named_gene2gene_dic.get(gene_id, gene_id)
    return taxa_dic.get(orig, "")


def get_leaf_voucher(leaf_name: str) -> str:
    """Extract voucher code from a leaf name.

    Args:
        leaf_name (str): Leaf name encoded as ``clade_voucher_index`` or ``voucher_index``.

    Returns:
        str: Voucher code parsed from the leaf name.

    Assumptions:
        Voucher codes are the second-to-last token when a clade prefix exists,
        or the first token when only ``voucher_index`` is present.
    """
    parts = leaf_name.split("_")
    if len(parts) >= 3:
        return parts[-2]
    if len(parts) == 2:
        return parts[0]
    return leaf_name


def build_leaf_annotations(
    Phylo_t: object,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
) -> tuple:
    """Annotate leaves with clade labels and voucher codes.

    Args:
        Phylo_t (object): ETE gene tree to annotate.
        taxa_dic (dict): Mapping from original gene IDs to clade labels.
        new_named_gene2gene_dic (dict): Mapping from renamed identifiers to original gene names.

    Returns:
        tuple: (leaf_to_clade, leaf_to_voucher) dictionaries.

    Assumptions:
        Leaf naming follows the ``clade_voucher_index`` convention.
    """
    leaf_to_clade = {}
    leaf_to_voucher = {}
    for leaf in Phylo_t.iter_leaves():
        leaf_to_clade[leaf.name] = get_leaf_clade(
            leaf.name,
            taxa_dic,
            new_named_gene2gene_dic,
        )
        leaf_to_voucher[leaf.name] = get_leaf_voucher(leaf.name)
    return leaf_to_clade, leaf_to_voucher


def compute_gene_tree_depths(Phylo_t: object) -> dict:
    """Compute topological depths for gene-tree nodes.

    Args:
        Phylo_t (object): ETE gene tree to traverse.

    Returns:
        dict: Mapping from node objects to integer depths from the root.

    Assumptions:
        The tree is connected and acyclic.
    """
    depths = {}
    for node in Phylo_t.traverse("preorder"):
        if node.is_root():
            depths[node] = 0
        else:
            depths[node] = depths[node.up] + 1
    return depths


def compute_species_tree_depths(species_tree: object) -> dict:
    """Compute topological depths for species-tree nodes.

    Args:
        species_tree (object): Rooted species tree.

    Returns:
        dict: Mapping from species-tree nodes to integer depths.

    Assumptions:
        The species tree is rooted and connected.
    """
    depths = {}
    for node in species_tree.traverse("preorder"):
        if node.is_root():
            depths[node] = 0
        else:
            depths[node] = depths[node.up] + 1
    return depths


def get_species_tree_depth_for_set(
    species_tree: object,
    depth_map: dict,
    species_set: set,
    cache: dict,
) -> int:
    """Return species-tree depth for the MRCA of a species set.

    Args:
        species_tree (object): Rooted species tree.
        depth_map (dict): Node-to-depth mapping for the species tree.
        species_set (set): Species identifiers mapped to voucher codes.
        cache (dict): Cache for previously computed species sets.

    Returns:
        int: Depth of the MRCA node, or 0 when mapping fails.

    Assumptions:
        Species identifiers exist as leaf names in the species tree.
    """
    if not species_set:
        return 0
    key = frozenset(species_set)
    if key in cache:
        return cache[key]
    try:
        if len(species_set) == 1:
            mapped_node = species_tree & list(species_set)[0]
        else:
            mapped_node = species_tree.get_common_ancestor(list(species_set))
        depth = depth_map.get(mapped_node, 0)
    except Exception:
        depth = 0
    cache[key] = depth
    return depth


def count_clade_leaves(node: object, leaf_to_clade: dict, target_clade: str) -> tuple:
    """Count target-clade and total leaves under a node.

    Args:
        node (object): ETE node to summarize.
        leaf_to_clade (dict): Mapping from leaf names to clade labels.
        target_clade (str): Target clade label.

    Returns:
        tuple: (target_count, total_count).

    Assumptions:
        Leaf annotations are precomputed and complete.
    """
    target_count = 0
    total_count = 0
    for leaf in node.iter_leaves():
        total_count += 1
        if leaf_to_clade.get(leaf.name) == target_clade:
            target_count += 1
    return target_count, total_count


def get_sister_species_for_target_clade(
    species_tree: object,
    leaf_to_clade: dict,
    leaf_to_voucher: dict,
    target_clade: str,
) -> set:
    """Get voucher species located on the immediate sister branch of a target clade.

    Args:
        species_tree (object): Rooted species tree.
        leaf_to_clade (dict): Mapping from gene-tree leaf names to clade labels.
        leaf_to_voucher (dict): Mapping from gene-tree leaf names to voucher species.
        target_clade (str): Target clade label currently being pruned.

    Returns:
        set: Voucher species set from the sister branch of the mapped target MRCA.

    Assumptions:
        The target clade is represented by at least one voucher in the current tree.
        If mapping fails or the target MRCA is the root, an empty set is returned.
    """
    target_species = {
        leaf_to_voucher[leaf_name]
        for leaf_name, clade_label in leaf_to_clade.items()
        if clade_label == target_clade and leaf_name in leaf_to_voucher
    }
    if not target_species:
        return set()

    try:
        if len(target_species) == 1:
            mapped_node = species_tree & list(target_species)[0]
        else:
            mapped_node = species_tree.get_common_ancestor(list(target_species))
    except Exception:
        return set()

    sister_species = set()
    for sister in mapped_node.get_sisters():
        sister_species.update(sister.get_leaf_names())
    return sister_species


def collect_dominant_lineages(
    Phylo_t: object,
    leaf_to_clade: dict,
    target_clade: str,
    dominant_purity: float,
) -> list:
    """Identify maximal dominant lineages for a target clade.

    Args:
        Phylo_t (object): ETE gene tree to scan.
        leaf_to_clade (dict): Mapping from leaf names to clade labels.
        target_clade (str): Target clade label.
        dominant_purity (float): Minimum purity threshold for dominant lineages.

    Returns:
        list: List of dominant lineage root nodes.

    Assumptions:
        A dominant lineage has more target tips than alien tips and purity above
        the specified cutoff, with at least two target tips.
    """
    dominant_roots = []
    processed_roots = set()
    # Preorder ensures we process ancestor candidates first; if any ancestor
    # has already been accepted as a dominant root, skip this node.
    for node in Phylo_t.traverse("preorder"):
        if any(ancestor in processed_roots for ancestor in node.get_ancestors()):
            continue
        target_count, total_count = count_clade_leaves(node, leaf_to_clade, target_clade)
        if target_count <= 1 or total_count == 0:
            continue
        alien_count = total_count - target_count
        purity = target_count / total_count
        if purity < dominant_purity or target_count <= alien_count:
            continue
        dominant_roots.append(node)
        processed_roots.add(node)
    return dominant_roots


def collect_alien_lineage_candidates(
    dominant_root: object,
    leaf_to_clade: dict,
    target_clade: str,
) -> list:
    """Collect alien lineage candidates within a dominant lineage.

    Args:
        dominant_root (object): Root node of a dominant lineage.
        leaf_to_clade (dict): Mapping from leaf names to clade labels.
        target_clade (str): Target clade label.

    Returns:
        list: Candidate nodes representing alien lineages.

    Assumptions:
        Alien families that are monophyletic are summarized by their MRCA;
        otherwise, individual tips are treated as separate candidates.
    """
    clade_to_leaves = {}
    for leaf in dominant_root.iter_leaves():
        clade_label = leaf_to_clade.get(leaf.name)
        if clade_label == target_clade:
            continue
        clade_to_leaves.setdefault(clade_label, []).append(leaf)

    candidates = []
    for clade_label, leaves in clade_to_leaves.items():
        if len(leaves) == 1:
            candidates.append(leaves[0])
            continue

        mrca = dominant_root.get_common_ancestor(leaves)
        if all(leaf_to_clade.get(name) == clade_label for name in mrca.get_leaf_names()):
            candidates.append(mrca)
        else:
            candidates.extend(leaves)

    return candidates


def normalize_values(values: list) -> list:
    """Normalize numeric values to the [0, 1] range.

    Args:
        values (list): Numeric values to normalize.

    Returns:
        list: Normalized values in the [0, 1] interval.

    Assumptions:
        When all values are identical, the normalized list is all zeros.
    """
    if not values:
        return []
    min_val = min(values)
    max_val = max(values)
    if max_val == min_val:
        return [0.0 for _ in values]
    return [(val - min_val) / (max_val - min_val) for val in values]


def score_alien_candidates(
    candidates: list,
    dominant_root: object,
    leaf_to_clade: dict,
    leaf_to_voucher: dict,
    species_tree: object,
    species_depths: dict,
    depth_cache: dict,
    gene_depths: dict,
    target_clade: str,
    sister_species_set: set,
) -> list:
    """Score alien lineage candidates within a dominant lineage.

    Args:
        candidates (list): Alien lineage candidate nodes.
        dominant_root (object): Root of the dominant lineage.
        leaf_to_clade (dict): Mapping from leaf names to clade labels.
        leaf_to_voucher (dict): Mapping from leaf names to voucher codes.
        species_tree (object): Species tree for distance estimation.
        species_depths (dict): Precomputed depth map for species-tree nodes.
        depth_cache (dict): Cache for MRCA depth lookups.
        gene_depths (dict): Precomputed depth map for gene-tree nodes.
        target_clade (str): Target clade label (e.g., Brassicaceae).
        sister_species_set (set): Voucher species on the immediate sister branch
            of the target clade in the species tree.

    Returns:
        list: List of candidate score dictionaries.

    Assumptions:
        Candidate lineages contain only alien tips.
    """
    dominant_leaves = dominant_root.get_leaf_names()
    dominant_size = len(dominant_leaves)
    brass_species_set = {
        leaf_to_voucher[name]
        for name in dominant_leaves
        if leaf_to_clade.get(name) == target_clade
    }
    brass_depth = get_species_tree_depth_for_set(
        species_tree,
        species_depths,
        brass_species_set,
        depth_cache,
    )

    results = []
    for node in candidates:
        leaf_names = node.get_leaf_names() if not node.is_leaf() else [node.name]
        alien_species_set = {leaf_to_voucher[name] for name in leaf_names}
        union_depth = get_species_tree_depth_for_set(
            species_tree,
            species_depths,
            brass_species_set | alien_species_set,
            depth_cache,
        )

        phylo_distance = brass_depth - union_depth
        alien_coverage = len(leaf_names) / dominant_size if dominant_size else 0.0
        alien_deepvar = gene_depths.get(node, 0) - gene_depths.get(dominant_root, 0)
        is_protected_sister = bool(
            alien_species_set and sister_species_set and alien_species_set.issubset(sister_species_set)
        )

        results.append(
            {
                "node": node,
                "leaf_names": leaf_names,
                "phylo_distance": phylo_distance,
                "alien_coverage": alien_coverage,
                "alien_deepvar": alien_deepvar,
                "is_protected_sister": is_protected_sister,
            }
        )

    phylo_values = [item["phylo_distance"] for item in results]
    deep_values = [item["alien_deepvar"] for item in results]
    phylo_norm = normalize_values(phylo_values)
    deep_norm = normalize_values(deep_values)

    for idx, item in enumerate(results):
        coverage_term = -np.log10(item["alien_coverage"] + 1e-4)
        item["combined_score"] = phylo_norm[idx] * deep_norm[idx] * coverage_term

    return results


def select_alien_tips_for_removal(
    scores: list,
    dominant_root: object,
    leaf_to_clade: dict,
    target_clade: str,
    final_purity: float,
    max_remove_fraction: float,
) -> set:
    """Select alien tips for removal based on ranked combined scores.

    Args:
        scores (list): Scored candidate dictionaries.
        dominant_root (object): Root of the dominant lineage.
        leaf_to_clade (dict): Mapping from leaf names to clade labels.
        target_clade (str): Target clade label.
        final_purity (float): Target purity for the dominant lineage after pruning.
        max_remove_fraction (float): Maximum fraction of tips removable from the lineage.

    Returns:
        set: Leaf names selected for removal.

    Assumptions:
        Dominant lineage purity is computed using the original lineage size.
    """
    target_count, total_count = count_clade_leaves(dominant_root, leaf_to_clade, target_clade)
    alien_count = total_count - target_count
    max_remove = max(int(max_remove_fraction * total_count), 1)
    alien_clade_counts = Counter(
        leaf_to_clade.get(leaf.name, "")
        for leaf in dominant_root.iter_leaves()
        if leaf_to_clade.get(leaf.name) != target_clade
    )

    removal_set = set()
    sorted_scores = sorted(scores, key=lambda x: x["combined_score"], reverse=True)

    for item in sorted_scores:
        if item.get("is_protected_sister", False):
            continue

        remaining_alien = max(alien_count - len(removal_set), 0)
        current_total = target_count + remaining_alien
        current_purity = target_count / current_total if current_total else 0.0

        if current_purity >= final_purity:
            break
        if len(removal_set) >= max_remove:
            break

        candidate_leaves = set(item["leaf_names"])
        new_leaves = candidate_leaves - removal_set
        if not new_leaves:
            continue
        # Keep at least one tip for every non-target clade in each dominant lineage.
        removed_by_clade = Counter(leaf_to_clade.get(name, "") for name in new_leaves)
        if any(
            clade_label
            and alien_clade_counts.get(clade_label, 0) <= removed_count
            for clade_label, removed_count in removed_by_clade.items()
        ):
            continue
        if len(removal_set) + len(new_leaves) > max_remove:
            break
        removal_set.update(new_leaves)

    return removal_set


# --------------------------
# 3. Pruning Core (Dominant-Lineage Based)
# --------------------------


def prune_all_clades(
    Phylo_t: object,
    species_tree: object,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
    final_purity: float,
    max_remove_fraction: float,
    log_handle,
    tre_ID: str,
    dominant_purity: float = 0.9,
    min_tips_per_clade: int = 1,
) -> object:
    t = Phylo_t.copy()
    leaf_to_clade, leaf_to_voucher = build_leaf_annotations(
        t, taxa_dic, new_named_gene2gene_dic
    )

    # Process only clades that are present in this tree and meet minimum size.
    clade_counts = Counter(leaf_to_clade.values())
    present_clades = [c for c, n in clade_counts.items() if n >= min_tips_per_clade]

    species_depths = compute_species_tree_depths(species_tree)
    depth_cache = {}
    gene_depths = compute_gene_tree_depths(t)

    global_remove = set()

    for target_clade in present_clades:
        remove_set = prune_one_clade_in_tree(
            t=t,
            species_tree=species_tree,
            leaf_to_clade=leaf_to_clade,
            leaf_to_voucher=leaf_to_voucher,
            new_named_gene2gene_dic=new_named_gene2gene_dic,
            species_depths=species_depths,
            depth_cache=depth_cache,
            gene_depths=gene_depths,
            target_clade=target_clade,
            dominant_purity=dominant_purity,
            final_purity=final_purity,
            max_remove_fraction=max_remove_fraction,
            log_handle=log_handle,
            tre_ID=tre_ID,
        )
        global_remove.update(remove_set)

    if global_remove:
        keep = set(t.get_leaf_names()) - global_remove
        t.prune(keep, preserve_branch_length=True)

    return t

def prune_one_clade_in_tree(
    t: object,
    species_tree: object,
    leaf_to_clade: dict,
    leaf_to_voucher: dict,
    new_named_gene2gene_dic: dict,
    species_depths: dict,
    depth_cache: dict,
    gene_depths: dict,
    target_clade: str,
    dominant_purity: float,
    final_purity: float,
    max_remove_fraction: float,
    log_handle,
    tre_ID: str,
) -> set:
    def to_original_gene_name(leaf_name: str) -> str:
        clade, gene_id = parse_leaf_components(leaf_name)
        if clade:
            return new_named_gene2gene_dic.get(gene_id, gene_id)
        parts = leaf_name.split("_")
        if len(parts) > 1:
            fallback_gene_id = "_".join(parts[1:])
            return new_named_gene2gene_dic.get(fallback_gene_id, fallback_gene_id)
        return new_named_gene2gene_dic.get(leaf_name, leaf_name)

    preorder_id = {node: idx for idx, node in enumerate(t.traverse("preorder"), start=1)}
    def get_tree_node_id(node: object) -> str:
        # Prefer the in-tree node numbering ID (e.g., N37) assigned by num_tre_node.
        return node.name if getattr(node, "name", "") else f"NODE_{preorder_id[node]}"

    dominant_roots = collect_dominant_lineages(
        t, leaf_to_clade, target_clade, dominant_purity
    )
    sister_species_set = get_sister_species_for_target_clade(
        species_tree=species_tree,
        leaf_to_clade=leaf_to_clade,
        leaf_to_voucher=leaf_to_voucher,
        target_clade=target_clade,
    )
    
    removal_set = set()
    for dominant_root in dominant_roots:
        candidates = collect_alien_lineage_candidates(
            dominant_root, leaf_to_clade, target_clade
        )
        if not candidates:
            continue

        scores = score_alien_candidates(
            candidates,
            dominant_root,
            leaf_to_clade,
            leaf_to_voucher,
            species_tree,
            species_depths,
            depth_cache,
            gene_depths,
            target_clade,
            sister_species_set,
        )

        selected = select_alien_tips_for_removal(
            scores,
            dominant_root,
            leaf_to_clade,
            target_clade,
            final_purity,
            max_remove_fraction,
        )

        dominant_target_count, dominant_total_count = count_clade_leaves(
            dominant_root, leaf_to_clade, target_clade
        )
        dominant_purity_val = (
            dominant_target_count / dominant_total_count if dominant_total_count else 0.0
        )
        dominant_root_leaf_count = len(dominant_root.get_leaf_names())

        # Record deletion rows only (one row per removed tip) without changing pruning behavior.
        for item in scores:
            candidate_node = item["node"]
            candidate_leaf_names = item["leaf_names"]

            candidate_leaf_set = set(candidate_leaf_names)
            removed_leaf_names = sorted(candidate_leaf_set & selected)
            if not removed_leaf_names:
                continue

            for removed_leaf in removed_leaf_names:
                removed_gene = to_original_gene_name(removed_leaf)
                log_handle.write(
                    f"{tre_ID}\t{target_clade}\t{get_tree_node_id(dominant_root)}\t"
                    f"{dominant_root_leaf_count}\t{dominant_target_count}\t{dominant_total_count}\t"
                    f"{dominant_purity_val}\t{removed_gene}\ttip\t1\t"
                    f"{item['phylo_distance']}\t{item['alien_coverage']}\t{item['alien_deepvar']}\t"
                    f"{item['combined_score']}\t1\t{removed_gene}\n"
                )

        removal_set.update(selected)

    return removal_set

# --------------------------
# 4. Visualization Utilities (Optional)
# --------------------------


def get_color_dict(taxa_dic: dict) -> dict:
    """Assign colors to taxa for visualization.

    Args:
        taxa_dic (dict): Mapping from gene identifiers to taxa labels.

    Returns:
        dict: Mapping from identifiers to ``taxon*color`` strings.

    Assumptions:
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
    """Apply consistent node styles and colored labels for visualization.

    Args:
        Phylo_t (object): ETE tree object to be styled.
        color_dict (dict): Mapping from gene identifiers to color strings.
        new_named_gene2gene_dic (dict): Mapping from renamed node identifiers to original gene names.

    Returns:
        object: The styled tree object.

    Assumptions:
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
            new_str = '_'.join(parts[1:]) if len(parts) > 1 else node.name
            gene = new_named_gene2gene_dic.get(new_str, new_str)
            color_info = color_dict.get(gene, None)
            if color_info:
                color = color_info.split("*")[-1]
                face = TextFace(species_name + "_" + gene, fgcolor=color, fstyle="italic")
                node.add_face(face, column=0)
    return Phylo_t


def generate_pdf(tre_ID: str, Phylo_t: object, tag: str) -> None:
    """Render a tree snapshot as a PDF file.

    Args:
        tre_ID (str): Tree identifier for output naming.
        Phylo_t (object): Tree to render.
        tag (str): Tag appended to the output filename.

    Returns:
        None

    Assumptions:
        ETE rendering is available in the runtime environment.
    """
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(f"{tre_ID}_{tag}", fsize=10), 0)
    Phylo_t.render(f"{tre_ID}_{tag}.pdf", tree_style=ts)


def merge_pdfs_side_by_side(f1: str, f2: str, out: str) -> None:
    """Merge two single-page PDFs horizontally.

    Args:
        f1 (str): Path to the first PDF.
        f2 (str): Path to the second PDF.
        out (str): Path to the output merged PDF.

    Returns:
        None

    Assumptions:
        Input PDFs contain at least one page.
    """
    writer = PdfWriter()
    with open(f1, "rb") as a, open(f2, "rb") as b:
        p1 = PdfReader(a).pages[0]
        p2 = PdfReader(b).pages[0]
        page = writer.add_blank_page(
            float(p1.mediabox.width) + float(p2.mediabox.width),
            max(float(p1.mediabox.height), float(p2.mediabox.height)),
        )
        page.merge_translated_page(p1, 0, 0)
        page.merge_translated_page(p2, float(p1.mediabox.width), 0)
        with open(out, "wb") as o:
            writer.write(o)


# --------------------------
# 5. Main Pipeline (Orchestration)
# --------------------------


def prune_main_Mono(
    tre_dic: dict,
    taxa_dic: dict,
    species_tree: object,
    purity_cutoff: float,
    max_remove_fraction: float,
    new_named_gene2gene_dic: dict,
    gene2new_named_gene_dic: dict,
    visual: bool = False,
) -> None:
    """Run the dominant-lineage mono-copy pruning pipeline.

    Args:
        tre_dic (dict): Mapping from tree IDs to file paths.
        taxa_dic (dict): Mapping from gene identifiers to clade labels.
        species_tree (object): Rooted species tree with voucher labels.
        purity_cutoff (float): Target purity for dominant lineages after pruning.
        max_remove_fraction (float): Maximum fraction of tips removable per lineage.
        new_named_gene2gene_dic (dict): Mapping from renamed identifiers to original gene names.
        gene2new_named_gene_dic (dict): Mapping from original gene names to renamed identifiers.
        visual (bool, optional): Whether to generate before/after PDF visualizations.

    Returns:
        None

    Assumptions:
        Input trees are valid and compatible with identifier mappings.
    """
    color_dic = get_color_dict(taxa_dic)

    base_dir = os.getcwd()
    module_root = (
        base_dir
        if os.path.basename(os.path.normpath(base_dir)) == "orthofilter_mono"
        else os.path.join(base_dir, "orthofilter_mono")
    )
    out_tree_dir = os.path.join(module_root, "pruned_tree")
    out_delete_gene_dir = os.path.join(module_root, "delete_gene")
    out_visual_dir = os.path.join(module_root, "visual")
    # Backward-compat cleanup: remove legacy duplicated insert logs.
    shutil.rmtree(os.path.join(module_root, "insert_gene"), ignore_errors=True)
    if visual:
        os.makedirs(out_visual_dir, exist_ok=True)
    else:
        out_visual_dir = None

    for d in (out_tree_dir, out_delete_gene_dir, out_visual_dir):
        if d is not None:
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d, exist_ok=True)


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
        with open(os.path.join(out_delete_gene_dir, f"{tre_ID}.delete_gene.tsv"), "w") as log:
            log.write(
                "tre_ID\ttarget_clade\tdominant_root_id\tdominant_root_leaf_count\t"
                "dominant_target_count\tdominant_total_count\tdominant_purity\t"
                "candidate_id\tcandidate_type\tcandidate_leaf_count\t"
                "phylo_distance\talien_coverage\t"
                "alien_deepvar\tcombined_score\tremoved_flag\tremoved_tips_original\n"
            )

            if visual:
                set_style(t, color_dic, new_named_gene2gene_dic)
                generate_pdf(tre_ID, t, "before")

            t1 = prune_all_clades(
                Phylo_t=t,
                species_tree=species_tree,
                taxa_dic=taxa_dic,
                new_named_gene2gene_dic=new_named_gene2gene_dic,
                final_purity=purity_cutoff,
                max_remove_fraction=max_remove_fraction,
                log_handle=log,
                tre_ID=tre_ID,
                dominant_purity=0.9,
                min_tips_per_clade=2,
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
        pbar.update(1)
    pbar.close()


# --------------------------
# 6. CLI Entry Point
# --------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Mono-copy ortholog pruning filter")
    parser.add_argument("--input_GF_list", required=True, help="Gene family list file")
    parser.add_argument("--input_taxa", required=True, help="Taxa mapping file")
    parser.add_argument("--input_sps_tree", required=True, help="Species tree file")
    parser.add_argument("--input_imap", default=None, help="Imap file for gene ID transfer (optional)")
    parser.add_argument("--purity_cutoff", type=float, default=0.95, help="Purity cutoff for clade assignment")
    parser.add_argument("--max_remove_fraction", type=float, default=0.5, help="Maximum fraction of tips to remove")
    parser.add_argument("--visual", action="store_true", help="Enable visual output")
    args = parser.parse_args()

    taxa_dic = read_and_return_dict(args.input_taxa)
    tre_dic = read_and_return_dict(args.input_GF_list)
    species_tree = PhyloTree(args.input_sps_tree)

    if args.input_imap:
        gene2new_named_gene_dic, new_named_gene2gene_dic, _, _ = gene_id_transfer(args.input_imap)
    else:
        gene2new_named_gene_dic = {}
        new_named_gene2gene_dic = {}

    prune_main_Mono(
        tre_dic,
        taxa_dic,
        species_tree=species_tree,
        purity_cutoff=args.purity_cutoff,
        max_remove_fraction=args.max_remove_fraction,
        new_named_gene2gene_dic=new_named_gene2gene_dic,
        gene2new_named_gene_dic=gene2new_named_gene_dic,
        visual=args.visual,
    )
