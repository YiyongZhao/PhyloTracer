"""
Gene duplication detection and reporting for the PhyloTracer pipeline.

This module identifies duplication clades in gene trees relative to a species
framework, classifies duplication patterns, and writes detailed summaries for
phylogenomic interpretation.
"""
from __future__ import annotations

import logging
from collections import Counter

logger = logging.getLogger(__name__)

import pandas as pd
from ete3 import PhyloTree
from tqdm import tqdm

from phylotracer import (
    read_phylo_tree,
    rename_input_tre,
    num_tre_node,
    annotate_gene_tree,
    find_dup_node,
    get_species_set,
    get_gene_pairs,
    gene_id_transfer,
    read_and_return_dict,
)

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
    gdtype_mode: str = "relaxed",
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
    logger.info("=== Gene Duplication Detection Configuration ===")
    logger.info("GD node support threshold: %s", duplication_support_threshold)
    logger.info("Subclade support threshold: %s", subclade_support_threshold)
    logger.info("Minimum duplicated species count: %s", duplicated_species_count_threshold)
    logger.info("Duplicated species percentage threshold: %s", duplicated_species_percentage_threshold)
    logger.info("Maximum variance of deepth: %s", max_topology_distance)
    logger.info("GD type assignment mode: %s", gdtype_mode)
    logger.info("===============================================")

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
                logger.warning("[Skip] %s is not a binary tree", tree_id)
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
            for clade in dup_node_list:
                species_set = get_species_set(species_tree&clade.map)
                child_a, child_b = clade.get_children()
                dup_species_count = len(get_species_set(child_a) & get_species_set(child_b))
                dup_ratio = dup_species_count / len(species_set) if species_set else 0
                gd_num_dict.setdefault(
                    voucher_to_taxa.get(clade.map, clade.map),
                    set(),
                ).add(clade)

                mapped_parent = species_tree & clade.map
                
                if gdtype_mode == "strict":
                    raw_model = get_model_strict(
                        clade,
                        species_tree,
                        max_topology_distance,
                    )
                else:
                    raw_model = get_model(clade, species_tree)
                gd_type_for_output = normalize_model(raw_model)
                
                if mapped_parent.is_leaf():
                    continue
                

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
    df.to_csv(f"gd_type_{gdtype_mode}.tsv", sep="\t", index=False)


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


def get_model_strict(clade: object, species_tree: object, deepvar: int) -> str:
    """
    Assign a strict raw GD model using branch-depth constrained matching.

    Parameters
    ----------
    clade : object
        Duplication node in the gene tree.
    species_tree : object
        Numbered species tree used for mapping.
    deepvar : int
        Maximum allowed absolute depth difference for A/B assignment.

    Returns
    -------
    str
        Raw four-character model string in ``A A B B`` slot order.

    Assumptions
    -----------
    The duplication node is binary and both child clades are binary.
    """
    map_node = species_tree & clade.map
    map_node_a, map_node_b = order_children_by_name(map_node)


    left_child, right_child = order_children_by_name(clade)
    if len(left_child.get_children()) != 2 or len(right_child.get_children()) != 2:
        return "XXXX"

    left_a, left_b = order_children_by_name(left_child)
    right_a, right_b = order_children_by_name(right_child)
    
    if abs(left_a.depth - map_node_a.depth) <= deepvar:
        type1='A'
    else:
        type1='X'
    
    if abs(right_a.depth - map_node_a.depth) <= deepvar:
        type2='A'
    else:
        type2='X'
    if abs(left_b.depth - map_node_b.depth) <= deepvar:
        type3='B'
    else:
        type3='X'
    
    if abs(right_b.depth - map_node_b.depth) <= deepvar:
        type4='B'
    else:
        type4='X'
    
    gdtype=type1+type2+type3+type4
    
    return gdtype


# ======================================================
# Section 3: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Gene duplication detection")
    parser.add_argument("--input_GF_list", required=True, help="Gene family list file")
    parser.add_argument("--input_imap", required=True, help="Imap file")
    parser.add_argument("--input_sps_tree", required=True, help="Species tree file")
    parser.add_argument("--gd_support", type=int, default=50, help="Duplication support threshold")
    parser.add_argument("--clade_support", type=int, default=50, help="Subclade support threshold")
    parser.add_argument("--dup_species_percent", type=float, default=0.5, help="Duplicated species percentage threshold")
    parser.add_argument("--dup_species_num", type=int, default=2, help="Duplicated species count threshold")
    parser.add_argument("--max_topology_distance", type=int, default=0, help="Maximum topology distance")
    parser.add_argument("--output", default="result.txt", help="Output file path")
    parser.add_argument("--gdtype_mode", default="relaxed", help="GD type mode (relaxed or strict)")
    args = parser.parse_args()

    gene_to_new_name, new_name_to_gene, voucher_to_taxa, _ = gene_id_transfer(args.input_imap)
    species_tree = PhyloTree(args.input_sps_tree)
    species_tree = rename_input_tre(species_tree, voucher_to_taxa)
    num_tre_node(species_tree)
    tree_paths = read_and_return_dict(args.input_GF_list)
    write_gene_duplication_results(
        args.output,
        tree_paths,
        args.gd_support,
        args.clade_support,
        args.dup_species_percent,
        args.dup_species_num,
        species_tree,
        gene_to_new_name,
        new_name_to_gene,
        voucher_to_taxa,
        args.max_topology_distance,
        gdtype_mode=args.gdtype_mode,
    )
