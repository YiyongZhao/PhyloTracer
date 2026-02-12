"""
Gene duplication detection and reporting for the PhyloTracer pipeline.

This module identifies duplication clades in gene trees relative to a species
framework, classifies duplication patterns, and writes detailed summaries for
phylogenomic interpretation.
"""

from calendar import c
from collections import Counter

from __init__ import *

# ======================================================
# Section 1: Duplication Detection and Reporting
# ======================================================


def write_gene_duplication_results(
    output_file: str,
    tree_paths: dict,
    duplication_support_threshold: int,
    subclade_support_threshold: int,
    duplicated_species_percentage_threshold: float,
    duplicated_species_count_threshold: int,
    species_tree: object,
    gene_to_new_name: dict,
    new_name_to_gene: dict,
    voucher_to_taxa: dict,
    max_topology_distance: int,
) -> None:
    """
    Detect gene duplications and write detailed and summary outputs.

    Parameters
    ----------
    output_file : str
        Path to the duplication result table.
    tree_paths : dict
        Mapping from tree identifiers to file paths.
    duplication_support_threshold : int
        Minimum support required for a duplication node.
    subclade_support_threshold : int
        Minimum support required for subclade validation.
    duplicated_species_percentage_threshold : float
        Minimum duplicated species proportion for GD classification.
    duplicated_species_count_threshold : int
        Minimum duplicated species count for GD classification.
    species_tree : object
        Species tree object used for reconciliation.
    gene_to_new_name : dict
        Mapping from original gene identifiers to renamed identifiers.
    new_name_to_gene : dict
        Mapping from renamed identifiers to original gene identifiers.
    voucher_to_taxa : dict
        Mapping from voucher identifiers to taxa labels.
    max_topology_distance : int
        Maximum allowed topological deviation for model classification.

    Returns
    -------
    None

    Assumptions
    -----------
    Gene trees are compatible with the species tree and ETE utilities.
    """
    print("=== Gene Duplication Detection Configuration ===")
    print(f"GD node support threshold: {duplication_support_threshold}")
    print(f"Subclade support threshold: {subclade_support_threshold}")
    print(f"Minimum duplicated species count: {duplicated_species_count_threshold}")
    print(
        "Duplicated species percentage threshold: "
        f"{duplicated_species_percentage_threshold}"
    )
    print(f"Maximum variance of deepth: {max_topology_distance}")
    print("===============================================")

    gd_type_dict: dict[str, list[str]] = {}
    gd_clade_count: dict[str, set[object]] = {}
    gd_num_dict: dict[str, set[object]] = {}


    gd_num = 1
    pbar = tqdm(total=len(tree_paths), desc="Processing trees", unit="tree")

    with open(output_file, "w") as fout:
        fout.write(
            "#tree_ID\tgd_id\tgd_support\tgene1\tgene2\t"
            "level\tspecies\tGD_dup_sps\tdup_ratio\tgd_type\tcomment\n"
        )

        for tree_id, tree_path in tree_paths.items():
            pbar.set_description(f"Processing {tree_id}")
            gene_tree = read_phylo_tree(tree_path)
            gene_tree = rename_input_tre(gene_tree, gene_to_new_name)
            
            if len(gene_tree.children) != 2:
                print(f"[Skip] {tree_id} is not a binary tree")
                continue

            num_tre_node(gene_tree)
            annotate_gene_tree(gene_tree, species_tree)
            
            dup_node_list = find_dup_node(
                gene_tree,
                species_tree,
                duplication_support_threshold,
                subclade_support_threshold,
                duplicated_species_count_threshold,
                duplicated_species_percentage_threshold,
                max_topology_distance,
            )


            for node in gene_tree.traverse("postorder"):
                if not hasattr(node, "map"):
                    continue
                parent = node.up
                if parent and hasattr(parent, "map") and parent.map == node.map:
                    continue
                level = voucher_to_taxa.get(node.map, node.map)
                if (species_tree & node.map).is_leaf():
                    continue

                gd_clade_count.setdefault(level, set()).add(node)
            vis_num=1
            for clade in dup_node_list:
                species_set = get_species_set(species_tree&clade.map)
                child_a, child_b = clade.get_children()
                dup_species_count = len(get_species_set(child_a) & get_species_set(child_b))
                dup_ratio = dup_species_count / len(species_set) if species_set else 0
                gd_num_dict.setdefault(
                    voucher_to_taxa.get(clade.map, clade.map),
                    set(),
                ).add(clade)

                
                if dup_species_count >=duplicated_species_count_threshold:
                    raw_model =get_model(clade, species_tree)
                if raw_model == 'Complex':
                    print(f"[Skip] {tree_id} {clade.map} is a complex model")
                mapped_parent = species_tree & clade.map
                
                if mapped_parent.is_leaf():
                    continue
                

                # overlap_sps = get_species_set(child_a) & get_species_set(child_b)
                # overlap_sps_mapped = map_species_set_to_node(species_tree, overlap_sps)



                # if clade.map=='S0':                    
                #     vis_clade=rename_input_tre(clade,new_name_to_gene)
                #     vis_clade.add_face(TextFace(clade.map, fsize=6, ftype="Arial",fgcolor='red'),column=0)
                #     vis_clade.add_face(TextFace(f"{clade.depth}", fsize=6, ftype="Arial",fgcolor='blue'),column=1,position="branch-bottom")
                #     a,b=vis_clade.get_children()
                #     a.add_face(TextFace(a.map, fsize=6, ftype="Arial",fgcolor='red'),column=0)
                #     a.add_face(TextFace(f"{a.depth}", fsize=6, ftype="Arial",fgcolor='blue'),column=1,position="branch-bottom")
                #     b.add_face(TextFace(b.map, fsize=6, ftype="Arial",fgcolor='red'),column=0)
                #     b.add_face(TextFace(f"{b.depth}", fsize=6, ftype="Arial",fgcolor='blue'),column=1,position="branch-bottom")
                #     a_a,a_b=a.get_children()
                #     b_a,b_b=b.get_children()

                #     a_a.add_face(TextFace(voucher_to_taxa.get(a_a.map, a_a.map), fsize=6, ftype="Arial",fgcolor='red'),column=0)
                #     a_a.add_face(TextFace(f"{a_a.depth}", fsize=6, ftype="Arial",fgcolor='blue'),column=1,position="branch-bottom")
                #     a_b.add_face(TextFace(f"{a_b.depth}", fsize=6, ftype="Arial",fgcolor='blue'),column=1,position="branch-bottom")
                #     a_b.add_face(TextFace(voucher_to_taxa.get(a_b.map, a_b.map), fsize=6, ftype="Arial",fgcolor='red'),column=0)
                #     b_a.add_face(TextFace(f"{b_a.depth}", fsize=6, ftype="Arial",fgcolor='blue'),column=1,position="branch-bottom")
                #     b_b.add_face(TextFace(f"{b_b.depth}", fsize=6, ftype="Arial",fgcolor='blue'),column=1,position="branch-bottom")
                #     b_a.add_face(TextFace(voucher_to_taxa.get(b_a.map, b_a.map), fsize=6, ftype="Arial",fgcolor='red'),column=0)
                #     b_b.add_face(TextFace(voucher_to_taxa.get(b_b.map, b_b.map), fsize=6, ftype="Arial",fgcolor='red'),column=0)
                    
                #     ts = TreeStyle()
                #     ts.scale = 10
                #     ts.legend_position = 1
                #     ts.show_leaf_name = False
                #     ts.guiding_lines_type = 0
                #     ts.guiding_lines_color = "black"
                #     ts.draw_guiding_lines = True
                #     ts.extra_branch_line_type = 0
                #     ts.extra_branch_line_color = "black"
                #     ts.legend.add_face(TextFace(f'Clade gd type is {raw_model}', fsize=6, ftype="Arial",fgcolor='red'), column=0)
                #     ts.legend.add_face(TextFace(f"Clade name is {clade.map} Overlap name is {overlap_sps_mapped.name}", fsize=6, ftype="Arial",fgcolor='red'), column=0)
                #     ts.legend.add_face(TextFace(f"Clade depth is {clade.depth} Overlap depth is {overlap_sps_mapped.depth}", fsize=6, ftype="Arial",fgcolor='red'), column=0)
                #     ts.legend.add_face(TextFace(f"Deepvar is ({abs(clade.depth-overlap_sps_mapped.depth)})", fsize=6, ftype="Arial",fgcolor='red'), column=0)

                #     vis_clade.render(f"{tree_id}_{vis_num}_{clade.map}.pdf", tree_style=ts)
                #     vis_num+=1

                



                gd_type_for_output = normalize_model(raw_model)
                
                gd_type_dict.setdefault(
                    voucher_to_taxa.get(clade.map, clade.map),
                    [],
                ).append(gd_type_for_output)

                gene_pairs = get_gene_pairs(clade)
                parent = clade.map if hasattr(clade, "map") else None
                level = voucher_to_taxa.get(parent, parent) if parent else "unknown"

                for sp, gene_a, gene_b in gene_pairs:
                    species = voucher_to_taxa.get(sp, sp)

                    fout.write(
                        f"{tree_id}\t{gd_num}\t{clade.support}\t"
                        f"{new_name_to_gene.get(gene_a, 'NA')}\t"
                        f"{new_name_to_gene.get(gene_b, 'NA')}\t"
                        f"{level}\t{species}\t"
                        f"{dup_species_count}\t{dup_ratio * 100:.2f}%\t"
                        f"{gd_type_for_output}\t-\n"
                    )

                gd_num += 1
            pbar.update()
        pbar.close()

    rows = []

    for node_name in gd_clade_count:
        gd_num = len(gd_num_dict.get(node_name, set()))
        clade_count = len(gd_clade_count[node_name])
        ratio = gd_num / clade_count if clade_count else 0
        ratio_str = f"{ratio * 100:.2f}%" if ratio else "0.00%"

        type_counter = Counter(gd_type_dict.get(node_name, []))

        rows.append(
            {
                "Newick_label": node_name,
                "GD": gd_num,
                "NUM": clade_count,
                "GD_ratio": ratio_str,
                "AABB": type_counter.get("AABB", 0),
                "AXBB": type_counter.get("AXBB", 0),
                "AABX": type_counter.get("AABX", 0),
                "Complex": type_counter.get("Complex", 0),
            }
        )

    df = pd.DataFrame(rows)
    df = df.sort_values("Newick_label", key=lambda x: x.str[1:].astype(int))
    df.to_csv("gd_type.tsv", sep="\t", index=False)


# ======================================================
# Section 2: Model Classification Utilities
# ======================================================


def normalize_model(raw_model: str) -> str:
    """
    Normalize a raw duplication model string into canonical types.

    Parameters
    ----------
    raw_model : str
        Raw model string composed of A/B/X labels.

    Returns
    -------
    str
        Canonical model label (ABAB, AAB, ABB, or Complex).

    Assumptions
    -----------
    Canonical models are defined by A/B counts and symmetry rules.
    """
    if raw_model in ("AXBB", "XABB"):
        return "AXBB"
    if raw_model in ("AABX", "AAXB"):
        return "AABX"
    if raw_model in ("AABB"):
        return "AABB"
    else:
        return 'Complex'

def order_children_by_name(n):
    c1, c2 = n.get_children()
    def key(x):
        m = x.num
        return int(m) if m else 10**18
    return (c1, c2) if key(c1) <= key(c2) else (c2, c1)

def get_model(clade: object, species_tree: object) -> str:
    """
    Assign GD type (ABAB / ABB / AAB / OTHER) in a tree2gd-equivalent manner.

    Notes
    -----
    - max_topology_distance is ignored (kept only for interface compatibility)
    - No topology distance is used
    - No grandchild decomposition is used
    """

    # 1. duplication node 在 species tree 上的映射节点 S
    map_node = species_tree & clade.map

    # 2. species tree 的二分：A / B
    map_node_A, map_node_B = order_children_by_name(map_node)
    A_species = get_species_set(map_node_A)
    B_species = get_species_set(map_node_B)

    # 3. duplication node 的两个子 clade（只看这一层）
    left_clade, right_clade = clade.get_children()

    # 4. 各子 clade 覆盖的物种集合
    left_species = get_species_set(left_clade)
    right_species = get_species_set(right_clade)

    gdtype=''
    gdtype+='A' if left_species & A_species else 'X'
    gdtype+='A' if right_species & A_species else 'X'
    gdtype+='B' if left_species & B_species else 'X'
    gdtype+='B' if right_species & B_species else 'X'

    return gdtype

# ======================================================
# Section 3: Gene Pair Extraction
# ======================================================


def get_gene_pairs(gd_node):
    """
    Extract gene pairs from a duplication node by species matching.

    Parameters
    ----------
    gd_node : object
        Duplication node whose child clades define gene pairs.

    Returns
    -------
    list
        List of (species, gene_left, gene_right) tuples.

    Assumptions
    -----------
    Leaf names are formatted as ``species_gene`` strings.
    """
    children = gd_node.get_children()
    if len(children) != 2:
        return []

    def collect_genes_by_species(clade):
        sp2genes = {}
        for leaf in clade.get_leaves():
            sp = leaf.name.split("_", 1)[0]
            if sp not in sp2genes:
                sp2genes[sp] = []
            sp2genes[sp].append(leaf.name)
        return sp2genes

    left_map = collect_genes_by_species(children[0])
    right_map = collect_genes_by_species(children[1])

    all_species = set(left_map.keys()) | set(right_map.keys())

    pairs = []
    for sp in all_species:
        left_genes = left_map.get(sp, [None])
        right_genes = right_map.get(sp, [None])

        for g1 in left_genes:
            for g2 in right_genes:
                pairs.append((sp, g1, g2))

    return pairs


# ======================================================
# Section 4: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    duplication_support_threshold = 50
    subclade_support_threshold = 50
    duplicated_species_percentage_threshold = 0.5
    duplicated_species_count_threshold = 2
    max_topology_distance = 0
    gene_to_new_name, new_name_to_gene, voucher_to_taxa = gene_id_transfer("imap.txt")
    species_tree = PhyloTree("30sptree.nwk")
    species_tree = rename_species_tree(species_tree, voucher_to_taxa)
    num_tre_node(species_tree)
    tree_paths = read_and_return_dict("GF.txt")
    output_file = "result.txt"
    write_gene_duplication_results(
        output_file,
        tree_paths,
        duplication_support_threshold,
        subclade_support_threshold,
        duplicated_species_percentage_threshold,
        duplicated_species_count_threshold,
        species_tree,
        gene_to_new_name,
        new_name_to_gene,
        voucher_to_taxa,
        max_topology_distance,
    )
