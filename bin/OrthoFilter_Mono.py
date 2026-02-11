"""
Monophyletic ortholog pruning guided by dominant-lineage purity in PhyloTracer.

This module detects dominant Brassicaceae lineages in gene trees and removes
alien tips using species-tree distance, coverage, and insertion depth metrics.
It optionally renders diagnostic PDFs and writes pruned trees for downstream
orthology inference.

Long-branch filtering is intentionally excluded and handled in a separate module.
"""

import matplotlib.colors as colors
import matplotlib.pyplot as plt
from PyPDF4 import PdfFileReader, PdfFileWriter

from __init__ import *
from BranchLength_NumericConverter import (
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
            node.name = f"{lineage}|{node.name}"
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
        parts = node.name.split("|")
        gene = parts[1]
        node.name = new_named_gene2gene_dic.get(gene, gene)
    return Phylo_t


# --------------------------
# 4. Dominant-Lineage Metrics
# --------------------------


def parse_leaf_components(leaf_name: str) -> tuple:
    """Parse clade, voucher, and gene index components from a leaf name.

    Args:
        leaf_name (str): Leaf name encoded as ``clade|voucher_index``.

    Returns:
        tuple: (clade_label, voucher_code, gene_index) with empty strings
            when components cannot be inferred.

    Assumptions:
        The final two vertical bar-delimited tokens represent the voucher and 
        gene index, while any preceding tokens belong to the clade label.
    """
    parts = leaf_name.split("|")
    if len(parts) >= 3:
        clade_label = "|".join(parts[:-2])
        voucher_code = parts[-2]
        gene_index = parts[-1]
        return clade_label, voucher_code, gene_index
    if len(parts) == 2:
        return parts[0], parts[1], ""
    return leaf_name, "", ""


def get_leaf_clade(
    leaf_name: str,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
) -> str:
    """Resolve the clade label for a leaf node.

    Args:
        leaf_name (str): Leaf name encoded as ``clade|voucher_index``.
        taxa_dic (dict): Mapping from original gene IDs to clade labels.
        new_named_gene2gene_dic (dict): Mapping from renamed identifiers to original gene names.

    Returns:
        str: Clade label for the leaf.

    Assumptions:
        Leaf names were prefixed with clade labels using ``rename_input_single_tre``.
    """
    clade_label, _, _ = parse_leaf_components(leaf_name)
    if clade_label and clade_label != leaf_name:
        return clade_label
    gene = leaf_name.split("|")[-1]
    return taxa_dic.get(gene, clade_label)


def get_leaf_voucher(leaf_name: str) -> str:
    """Extract voucher code from a leaf name.

    Args:
        leaf_name (str): Leaf name encoded as ``clade|voucher_index`` or ``voucher_index``.

    Returns:
        str: Voucher code parsed from the leaf name.

    Assumptions:
        Voucher codes are the second-to-last token when a clade prefix exists,
        or the first token when only ``voucher_index`` is present.
    """
    parts = leaf_name.split("|")
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
        Leaf naming follows the ``clade|voucher_index`` convention.
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
    for node in Phylo_t.traverse("preorder"):
        target_count, total_count = count_clade_leaves(node, leaf_to_clade, target_clade)
        if target_count <= 1 or total_count == 0:
            continue
        alien_count = total_count - target_count
        purity = target_count / total_count
        if purity < dominant_purity or target_count <= alien_count:
            continue
        if any(ancestor in dominant_roots for ancestor in node.get_ancestors()):
            continue
        dominant_roots.append(node)
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

        results.append(
            {
                "node": node,
                "leaf_names": leaf_names,
                "phylo_distance": phylo_distance,
                "alien_coverage": alien_coverage,
                "alien_deepvar": alien_deepvar,
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

    removal_set = set()
    sorted_scores = sorted(scores, key=lambda x: x["combined_score"], reverse=True)

    for item in sorted_scores:
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
        if len(removal_set) + len(new_leaves) > max_remove:
            break
        removal_set.update(new_leaves)

    return removal_set


# --------------------------
# 5. Pruning Core (Dominant-Lineage Based)
# --------------------------


def prune_brassicaceae_lineages(
    Phylo_t: object,
    species_tree: object,
    taxa_dic: dict,
    new_named_gene2gene_dic: dict,
    purity_cutoff: float,
    max_remove_fraction: float,
    log_handle,
    tre_ID: str,
) -> object:
    """Prune alien tips from Brassicaceae dominant lineages.

    Args:
        Phylo_t (object): ETE gene tree to prune.
        species_tree (object): Rooted species tree with voucher labels.
        taxa_dic (dict): Mapping from original gene IDs to clade labels.
        new_named_gene2gene_dic (dict): Mapping from renamed identifiers to original gene names.
        purity_cutoff (float): Target purity for dominant lineages after pruning.
        max_remove_fraction (float): Maximum fraction of tips removable per lineage.
        log_handle (object): File handle for logging candidate scores.
        tre_ID (str): Tree identifier for logging.

    Returns:
        object: Pruned gene tree.

    Assumptions:
        Brassicaceae tips are encoded via clade-prefixed leaf names.
    """
    target_clade = "Brassicaceae"
    dominant_purity = 0.9

    t = Phylo_t.copy()
    leaf_to_clade, leaf_to_voucher = build_leaf_annotations(
        t,
        taxa_dic,
        new_named_gene2gene_dic,
    )

    species_depths = compute_species_tree_depths(species_tree)
    depth_cache = {}
    gene_depths = compute_gene_tree_depths(t)

    dominant_roots = collect_dominant_lineages(
        t,
        leaf_to_clade,
        target_clade,
        dominant_purity,
    )

    removal_set = set()
    for dominant_root in dominant_roots:
        candidates = collect_alien_lineage_candidates(
            dominant_root,
            leaf_to_clade,
            target_clade,
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
        )

        selected = select_alien_tips_for_removal(
            scores,
            dominant_root,
            leaf_to_clade,
            target_clade,
            purity_cutoff,
            max_remove_fraction,
        )

        for item in scores:
            candidate_name = item["node"].name if item["node"].name else "NA"
            remove_flag = 1 if set(item["leaf_names"]).issubset(selected) else 0
            log_handle.write(
                f"{tre_ID}\t{dominant_root.name}\t{candidate_name}\t"
                f"{item['phylo_distance']}\t{item['alien_coverage']}\t"
                f"{item['alien_deepvar']}\t{item['combined_score']}\t{remove_flag}\n"
            )

        removal_set.update(selected)

    if removal_set:
        keep = set(t.get_leaf_names()) - removal_set
        t.prune(keep, preserve_branch_length=True)

    return t


# --------------------------
# 6. Visualization Utilities (Optional)
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
            parts = node.name.split("|")
            species_name = parts[0]
            new_str = parts[1] if len(parts) > 1 else node.name
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


# --------------------------
# 7. Main Pipeline (Orchestration)
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
        log.write(
            "tre_ID\tdominant_root\tcandidate\tphylo_distance\t"
            "alien_coverage\talien_deepvar\tcombined_score\tremoved\n"
        )

        if visual:
            set_style(t, color_dic, new_named_gene2gene_dic)
            generate_pdf(tre_ID, t, "before")

        t1 = t
        if len(get_species_set(t)) > 1:
            t1 = prune_brassicaceae_lineages(
                t,
                species_tree,
                taxa_dic,
                new_named_gene2gene_dic,
                purity_cutoff,
                max_remove_fraction,
                log,
                tre_ID,
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
    
        log.close()
    pbar.close()


# --------------------------
# 8. CLI Entry Point
# --------------------------

if __name__ == "__main__":
    taxa_dic = read_and_return_dict("taxa.txt")
    tre_dic = read_and_return_dict("100_nosingle_GF_list.txt")
    species_tree = PhyloTree("sps_tree.txt")
    prune_main_Mono(
        tre_dic,
        taxa_dic,
        species_tree=species_tree,
        purity_cutoff=0.95,
        max_remove_fraction=0.5,
        new_named_gene2gene_dic={},
        gene2new_named_gene_dic={},
    )
