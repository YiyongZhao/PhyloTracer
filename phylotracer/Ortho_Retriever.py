"""
Ortholog retrieval and single-copy tree extraction for PhyloTracer analyses.

This module identifies offcut paralogous clades, separates ortholog groups,
extracts single-copy subtrees, and writes summary tables for downstream use.
"""

import os
import re
from collections import defaultdict

from ete3 import PhyloTree

from phylotracer import (
    find_tre_dup,
    gene_id_transfer,
    is_rooted,
    read_and_return_dict,
    read_phylo_tree,
    rename_input_tre,
    root_tre_with_midpoint_outgroup,
    serialize_tree_by_input_branch_length_style,
)

# ======================================================
# Section 1: Identifier and Tree Preparation Utilities
# ======================================================


def rename_len_dic(len_dic: dict, gene2new_named_gene_dic: dict) -> dict:
    """
    Rename length dictionary keys using the renamed gene mapping.

    Parameters
    ----------
    len_dic : dict
        Mapping from original gene identifiers to sequence lengths.
    gene2new_named_gene_dic : dict
        Mapping from original gene identifiers to renamed identifiers.

    Returns
    -------
    dict
        Mapping from renamed identifiers to sequence lengths.

    Assumptions
    -----------
    Only genes present in ``gene2new_named_gene_dic`` are retained.
    """
    return {
        gene2new_named_gene_dic[key]: value
        for key, value in len_dic.items()
        if key in gene2new_named_gene_dic.keys()
    }


def count_sps_num(ev_seqs: set) -> int:
    """
    Count unique species represented by gene identifiers.

    Parameters
    ----------
    ev_seqs : set
        Set of gene identifiers.

    Returns
    -------
    int
        Number of unique species inferred from gene prefixes.

    Assumptions
    -----------
    Species codes are encoded in the first three characters of gene names.
    """
    sps_set = {gene[0:3] for gene in ev_seqs}
    return len(sps_set)


def extract_tree(ev_seqs: set, Phylo_t: object) -> object:
    """
    Prune a tree to retain only a specified gene set.

    Parameters
    ----------
    ev_seqs : set
        Set of gene identifiers to retain.
    Phylo_t : object
        ETE tree object to be pruned.

    Returns
    -------
    object
        Pruned tree containing only the specified genes.

    Assumptions
    -----------
    All identifiers in ``ev_seqs`` exist in the tree leaves.
    """
    Phylo_t.prune(ev_seqs)
    return Phylo_t


# ======================================================
# Section 2: Offcut Detection and Splitting
# ======================================================


def offcut_tre(Phylo_t: object, renamed_len_dic: dict) -> list:
    """
    Identify offcut gene sets based on duplication structure and lengths.

    Parameters
    ----------
    Phylo_t : object
        ETE tree object representing a gene family.
    renamed_len_dic : dict
        Mapping from renamed gene identifiers to sequence lengths.

    Returns
    -------
    list
        Tuple of (principal_gene_set, filtered_offcut_sets).

    Assumptions
    -----------
    Duplication structures are encoded in ``find_tre_dup`` output.
    """
    tre_ParaL, GF_leaves_S = find_tre_dup(Phylo_t)
    offcut_ev_seqs_L0 = []
    for ev_seqs in tre_ParaL:
        fd = ev_seqs.strip().split("<=>")
        if len(fd) != 2:  # guard against malformed input without <=> separator
            continue
        ev_seqs1, ev_seqs2 = set(fd[0].split(",")), set(fd[1].split(","))
        if None in ev_seqs1 or None in ev_seqs2:
            continue
        sps_num1 = count_sps_num(ev_seqs1)
        seq_num1 = len(ev_seqs1)
        sps_num2 = count_sps_num(ev_seqs2)
        seq_num2 = len(ev_seqs2)
        if sps_num1 == sps_num2:
            if seq_num1 > seq_num2:
                offcut_ev_seqs_L0.append(ev_seqs2)
            elif seq_num1 < seq_num2:
                offcut_ev_seqs_L0.append(ev_seqs1)
            elif seq_num1 == seq_num2:
                if renamed_len_dic:
                    ev_seqs1_avg_len = (
                        sum(
                            [
                                renamed_len_dic.get(Gene, 0)
                                for Gene in ev_seqs1
                                if Gene is not None
                            ]
                        )
                        / sps_num1
                        if sps_num1
                        else 0
                    )
                    ev_seqs2_avg_len = (
                        sum(
                            [
                                renamed_len_dic.get(Gene, 0)
                                for Gene in ev_seqs2
                                if Gene is not None
                            ]
                        )
                        / sps_num2
                        if sps_num2
                        else 0
                    )
                else:
                    ev_seqs1_avg_len = 0
                    ev_seqs2_avg_len = 0
                if ev_seqs1_avg_len >= ev_seqs2_avg_len:
                    offcut_ev_seqs_L0.append(ev_seqs2)
                else:
                    offcut_ev_seqs_L0.append(ev_seqs1)
        else:
            if sps_num1 > sps_num2:
                offcut_ev_seqs_L0.append(ev_seqs2)
            else:
                offcut_ev_seqs_L0.append(ev_seqs1)
    offcut_Gene_S = set()
    for ev_seqs in offcut_ev_seqs_L0:
        offcut_Gene_S = offcut_Gene_S | ev_seqs
    principal_gene_S = GF_leaves_S - offcut_Gene_S
    filtered_offcut_ev_seqs_L0 = []
    for ev_seqs in offcut_ev_seqs_L0:
        sps_num = count_sps_num(ev_seqs)
        if sps_num >= 2:
            filtered_offcut_ev_seqs_L0.append(ev_seqs)
    filtered_offcut_ev_seqs_L0 = rm_dup(filtered_offcut_ev_seqs_L0)
    return principal_gene_S, filtered_offcut_ev_seqs_L0


def rm_dup(paralogs_L):
    """
    Remove redundant paralog sets that are subsets of others.

    Parameters
    ----------
    paralogs_L : list
        List of paralogous gene sets.

    Returns
    -------
    list
        Deduplicated list of paralog sets.

    Assumptions
    -----------
    Gene sets are represented as Python sets.
    """
    result = []
    for i, ev_seqs1 in enumerate(paralogs_L):
        is_proper_subset = False
        for j, ev_seqs2 in enumerate(paralogs_L):
            if i != j and ev_seqs1 < ev_seqs2:
                is_proper_subset = True
                break
        if not is_proper_subset:
            result.append(ev_seqs1)
    return result


def split_offcut_ev_seqs(offcut_ev_seqs_L0: list) -> list:
    """
    Split offcut sets into ortholog and paralog collections.

    Parameters
    ----------
    offcut_ev_seqs_L0 : list
        List of offcut gene sets.

    Returns
    -------
    list
        Tuple of (ortholog_sets, paralog_sets).

    Assumptions
    -----------
    Species counts distinguish orthologs from paralogs.
    """
    othologs_L = []
    paralogs_L = []
    for ev_seqs in offcut_ev_seqs_L0:
        sps_num, seq_num = count_sps_num(ev_seqs), len(ev_seqs)
        if sps_num >= 2:
            if seq_num > sps_num:
                paralogs_L.append(ev_seqs)
            else:
                othologs_L.append(ev_seqs)
    paralogs_L = rm_dup(paralogs_L)
    return othologs_L, paralogs_L


def iterator(
    offcut_ev_seqs_L0: list,
    Phylo_t: object,
    new_named_gene2gene_dic: dict,
    minor_othologs_L: list,
    tre_path: str,
    renamed_len_dic: dict,
) -> list:
    """
    Recursively extract ortholog sets from offcut paralogs.

    Parameters
    ----------
    offcut_ev_seqs_L0 : list
        List of offcut gene sets.
    Phylo_t : object
        ETE tree object for the gene family.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene identifiers.
    minor_othologs_L : list
        Accumulator list of minor ortholog sets.
    tre_path : str
        Path to the gene tree file.
    renamed_len_dic : dict
        Mapping from renamed gene identifiers to sequence lengths.

    Returns
    -------
    list
        Updated list of minor ortholog sets.

    Assumptions
    -----------
    Recursive pruning terminates when no paralogs remain.
    """
    othologs_L, paralogs_L = split_offcut_ev_seqs(offcut_ev_seqs_L0)
    minor_othologs_L += othologs_L
    if paralogs_L != []:
        for i, ev_seqs in enumerate(paralogs_L):
            Phylo_t = read_phylo_tree(tre_path)
            if is_rooted(Phylo_t):
                pass
            else:
                Phylo_t = root_tre_with_midpoint_outgroup(Phylo_t)
            Phylo_t = rename_input_tre(Phylo_t, new_named_gene2gene_dic)
            Phylo_t = extract_tree(ev_seqs, Phylo_t)
            principal_gene_S, offcut_ev_seqs_L0 = offcut_tre(Phylo_t, renamed_len_dic)
            minor_othologs_L.append(principal_gene_S)
            iterator(
                offcut_ev_seqs_L0,
                Phylo_t,
                new_named_gene2gene_dic,
                minor_othologs_L,
                tre_path,
                renamed_len_dic,
            )
    return minor_othologs_L


# ======================================================
# Section 3: Ortholog Group Naming and Extraction
# ======================================================


def rename_OGs_tre_name(
    principal_gene_S: list,
    minor_othologs_L: list,
    tre_ID: str,
) -> list:
    """
    Assign ordered names to ortholog groups for output tracking.

    Parameters
    ----------
    principal_gene_S : list
        Main ortholog gene set.
    minor_othologs_L : list
        List of minor ortholog gene sets.
    tre_ID : str
        Tree identifier used as a prefix.

    Returns
    -------
    list
        List of (tree_name, gene_set) tuples.

    Assumptions
    -----------
    Ortholog group sizes determine output ordering.
    """
    OG_count = len(minor_othologs_L) + 1
    tre_name_L = ["T" + str(i + 1) for i in range(OG_count)]
    sps_num_L = [len(OG) for OG in minor_othologs_L]

    tre_name2sps_num_dict = dict(zip(tre_name_L, sps_num_L))
    tre_name2OG_dict = dict(zip(tre_name_L, minor_othologs_L))
    orderd_OG_sps_num_L = sorted(
        tre_name2sps_num_dict.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    ordered_name_OG_L = []
    ordered_name_OG_L.append(
        (str(tre_ID) + "_" + tre_name_L[0] + "_" + str(len(principal_gene_S)), principal_gene_S)
    )
    index = 0
    for name, count in orderd_OG_sps_num_L:
        index += 1
        ordered_name_OG_L.append(
            (str(tre_ID) + "_" + tre_name_L[index] + "_" + str(count), tre_name2OG_dict[name])
        )
    return ordered_name_OG_L


def get_single_copy_trees(
    Phylo_t1: object,
    renamed_len_dic: dict,
    gene2new_named_gene_dic: dict,
    new_named_gene2gene_dic: dict,
    tre_path: str,
    tre_ID: str,
) -> list:
    """
    Extract single-copy ortholog subtrees from a gene family tree.

    Parameters
    ----------
    Phylo_t1 : object
        Gene family tree object.
    renamed_len_dic : dict
        Mapping from renamed gene identifiers to sequence lengths.
    gene2new_named_gene_dic : dict
        Mapping from original gene identifiers to renamed identifiers.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene identifiers.
    tre_path : str
        Path to the gene tree file.
    tre_ID : str
        Tree identifier used in output naming.

    Returns
    -------
    list
        List of (tree_name, subtree_object) tuples.

    Assumptions
    -----------
    Input tree identifiers are consistent across mapping dictionaries.
    """
    trees = []
    principal_gene_S, filtered_offcut_ev_seqs_L0 = offcut_tre(Phylo_t1, renamed_len_dic)
    minor_othologs_L = []
    minor_othologs_L = iterator(
        filtered_offcut_ev_seqs_L0,
        Phylo_t1,
        gene2new_named_gene_dic,
        minor_othologs_L,
        tre_path,
        renamed_len_dic,
    )
    ordered_name_OG_L = rename_OGs_tre_name(principal_gene_S, minor_othologs_L, tre_ID)
    for tre_name, OG_S in ordered_name_OG_L:
        OG_L = [new_named_gene2gene_dic.get(OG, None) for OG in OG_S]
        if None in OG_L:
            continue
        Phylo_t0 = read_phylo_tree(tre_path)
        if is_rooted(Phylo_t0):
            Phylo_t = Phylo_t0
        else:
            Phylo_t = root_tre_with_midpoint_outgroup(Phylo_t0)
        Phylo_t_OG_L = extract_tree(OG_L, Phylo_t)
        trees.append((tre_name, Phylo_t_OG_L))
    return trees


def parse_synteny_blocks(synteny_blocks_path: str) -> dict:
    """
    Parse a raw synteny block file into a gene-to-block mapping.

    Parameters
    ----------
    synteny_blocks_path : str
        Path to a raw block file where ``#`` starts a new block and each
        non-comment line stores a gene pair.

    Returns
    -------
    dict
        Mapping from gene identifier to a set of block IDs.

    Assumptions
    -----------
    Each non-comment line contains either:
    1. a simplified two-column gene pair, or
    2. a WGDI-style row where the first and third columns are gene IDs.
    """
    gene_to_blocks = defaultdict(set)
    current_block = None
    block_index = 0

    with open(synteny_blocks_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                block_index += 1
                header = line.lstrip("#").strip()
                block_label = ""
                if header:
                    parts = header.split()
                    if len(parts) >= 2 and parts[0] == "Alignment":
                        block_label = f"{parts[0]}_{parts[1].rstrip(':')}"
                    else:
                        block_label = parts[0]
                current_block = block_label or f"block_{block_index}"
                continue

            if current_block is None:
                block_index += 1
                current_block = f"block_{block_index}"

            parts = line.split()
            if len(parts) < 2:
                continue
            gene1 = parts[0]
            gene2 = parts[1]

            # Support WGDI-like rows such as:
            # gene1  order1  gene2  order2  strand
            # while keeping the original two-column format unchanged.
            if len(parts) >= 3 and parts[1].lstrip("-").isdigit():
                gene2 = parts[2]

            gene_to_blocks[gene1].add(current_block)
            gene_to_blocks[gene2].add(current_block)

    return gene_to_blocks


def filter_trees_by_synteny(trees: list, gene_to_blocks: dict) -> tuple[list, list]:
    """
    Filter candidate ortholog trees by raw synteny block support.

    Parameters
    ----------
    trees : list
        List of ``(tree_name, tree_object)`` tuples.
    gene_to_blocks : dict
        Mapping from gene identifier to a set of synteny block IDs.

    Returns
    -------
    tuple[list, list]
        Filtered trees and per-tree synteny report rows.

    Assumptions
    -----------
    A candidate passes if at least two genes have synteny support and the most
    supported block covers at least half of the synteny-supported genes.
    """
    filtered_trees = []
    report_rows = []

    for tree_name, phylo_tree in trees:
        leaf_names = [leaf.name for leaf in phylo_tree]
        block_gene_counts = defaultdict(set)
        supported_genes = 0

        for gene in leaf_names:
            blocks = gene_to_blocks.get(gene, set())
            if not blocks:
                continue
            supported_genes += 1
            for block_id in blocks:
                block_gene_counts[block_id].add(gene)

        best_block = ""
        best_block_gene_count = 0
        if block_gene_counts:
            best_block, best_genes = max(
                block_gene_counts.items(),
                key=lambda item: (len(item[1]), item[0]),
            )
            best_block_gene_count = len(best_genes)

        support_ratio = (
            best_block_gene_count / supported_genes if supported_genes else 0.0
        )
        passes = supported_genes >= 2 and support_ratio >= 0.5
        if passes:
            filtered_trees.append((tree_name, phylo_tree))

        report_rows.append(
            (
                tree_name,
                str(len(leaf_names)),
                str(supported_genes),
                best_block,
                str(best_block_gene_count),
                f"{support_ratio:.4f}",
                "pass" if passes else "fail",
            )
        )

    return filtered_trees, report_rows


# ======================================================
# Section 4: Main Pipeline (Orchestration)
# ======================================================

TREE_NAME_PATTERN = re.compile(r"^(?P<gf_id>.+)_T\d+_\d+$")


def parse_original_tree_id(tree_name: str) -> str:
    """Extract the original gene-family ID from an Ortho_Retriever tree name."""
    match = TREE_NAME_PATTERN.match(tree_name.strip())
    if not match:
        raise ValueError(f"Cannot parse original tree ID from tree name: {tree_name}")
    return match.group("gf_id")


def choose_strict_single_copy_outgroup_genes(
    candidate_genes: list[str],
    gene_to_species: dict[str, str],
    ingroup_species: set[str],
) -> list[str]:
    """
    Keep only strict outgroup genes from species absent from the ingroup.

    A species is retained only when it contributes exactly one gene in the
    candidate sister clade. Multi-copy candidate species are rejected rather
    than collapsed to an arbitrary representative.
    """
    genes_by_species: dict[str, list[str]] = {}
    for gene in candidate_genes:
        species = gene_to_species.get(gene)
        if species is None or species in ingroup_species:
            continue
        genes_by_species.setdefault(species, []).append(gene)

    selected = []
    for species in sorted(genes_by_species):
        species_genes = genes_by_species[species]
        if len(species_genes) == 1:
            selected.append(species_genes[0])
    return sorted(selected)


def get_preferred_sister_subclade(mrca_node) -> object | None:
    """Return the largest child clade on the sister side of the ingroup MRCA."""
    parent = mrca_node.up
    if parent is None:
        return None

    sister_nodes = [child for child in parent.children if child is not mrca_node]
    if not sister_nodes:
        return None

    sister = max(
        sister_nodes,
        key=lambda node: (len(node.get_leaf_names()), node.name),
    )
    if sister.is_leaf():
        return sister

    child_clades = list(sister.children)
    if not child_clades:
        return sister

    return max(
        child_clades,
        key=lambda node: (len(node.get_leaf_names()), node.name),
    )


def select_outgroup_genes(
    ingroup_genes: list[str],
    original_tree,
    gene_to_species: dict[str, str],
) -> list[str]:
    """
    Select strict outgroup genes for one ortholog group.

    The outgroup must come from the sister side of the ingroup MRCA, must not
    share species with the ingroup, and each retained outgroup species must be
    single-copy within the candidate sister clade.
    """
    if not ingroup_genes:
        return []

    if len(ingroup_genes) == 1:
        mrca_node = original_tree & ingroup_genes[0]
    else:
        mrca_node = original_tree.get_common_ancestor(ingroup_genes)

    if mrca_node.is_root():
        return []

    sister_clade = get_preferred_sister_subclade(mrca_node)
    if sister_clade is None:
        return []

    if sister_clade.is_leaf():
        candidate_genes = [sister_clade.name]
    else:
        candidate_genes = sister_clade.get_leaf_names()

    ingroup_gene_set = set(ingroup_genes)
    candidate_genes = [gene for gene in candidate_genes if gene not in ingroup_gene_set]
    if not candidate_genes:
        return []

    ingroup_species = {
        gene_to_species[gene]
        for gene in ingroup_genes
        if gene in gene_to_species
    }
    return choose_strict_single_copy_outgroup_genes(
        candidate_genes,
        gene_to_species,
        ingroup_species,
    )


def build_rooted_combined_tree(original_tree, ingroup_genes: list[str], outgroup_genes: list[str]):
    """Prune the original tree to ingroup + outgroup genes and reroot by outgroup."""
    combined = original_tree.copy("newick")
    keep_genes = list(dict.fromkeys(ingroup_genes + outgroup_genes))
    combined.prune(keep_genes)

    if outgroup_genes:
        if len(outgroup_genes) == 1:
            outgroup_node = combined & outgroup_genes[0]
        else:
            outgroup_node = combined.get_common_ancestor(outgroup_genes)
        combined.set_outgroup(outgroup_node)
    return combined


def split_main(
    tre_dic,
    gene2new_named_gene_dic,
    new_named_gene2gene_dic,
    renamed_len_dic,
    synteny_blocks_path: str | None = None,
    add_outgroup: bool = False,
    gene_to_species: dict[str, str] | None = None,
):
    """
    Split gene family trees into single-copy ortholog trees and write outputs.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to file paths.
    gene2new_named_gene_dic : dict
        Mapping from original gene identifiers to renamed identifiers.
    new_named_gene2gene_dic : dict
        Mapping from renamed identifiers to original gene identifiers.
    renamed_len_dic : dict
        Mapping from renamed gene identifiers to sequence lengths.

    Returns
    -------
    None

    Assumptions
    -----------
    Input trees can be rooted and pruned without identifier conflicts.
    """
    import csv

    gene_to_blocks = {}
    if synteny_blocks_path:
        gene_to_blocks = parse_synteny_blocks(synteny_blocks_path)

    outgroup_report_rows = []

    with open("ortho_retriever_summary.txt", "w") as o:
        o.write("tre_name" + "\t" + "single_copy_tree" + "\n")

        processed_lines = []
        tsv_data = {}
        synteny_report_rows = []

        for tre_ID, tre_path in tre_dic.items():
            Phylo_t0 = read_phylo_tree(tre_path)
            if len(Phylo_t0) == 0 or len(Phylo_t0.children) == 0:
                continue
            if is_rooted(Phylo_t0):
                Phylo_t1 = Phylo_t0
            else:
                Phylo_t1 = root_tre_with_midpoint_outgroup(Phylo_t0)

            Phylo_t2 = rename_input_tre(Phylo_t1, gene2new_named_gene_dic)
            trees = get_single_copy_trees(
                Phylo_t2,
                renamed_len_dic,
                gene2new_named_gene_dic,
                new_named_gene2gene_dic,
                tre_path,
                tre_ID,
            )
            if synteny_blocks_path:
                trees, tree_report_rows = filter_trees_by_synteny(trees, gene_to_blocks)
                synteny_report_rows.extend((tre_ID,) + row for row in tree_report_rows)

            tree_ortholog_trees = []
            for clade in trees:
                tre_name = clade[0]
                Phylo_t_OG_L = clade[1]
                tree_to_write = Phylo_t_OG_L
                if add_outgroup and gene_to_species is not None:
                    original_tree = read_phylo_tree(tre_path)
                    if not is_rooted(original_tree):
                        original_tree = root_tre_with_midpoint_outgroup(original_tree)
                    ingroup_genes = [leaf.name for leaf in Phylo_t_OG_L]
                    outgroup_genes = select_outgroup_genes(
                        ingroup_genes,
                        original_tree,
                        gene_to_species,
                    )
                    if outgroup_genes:
                        tree_to_write = build_rooted_combined_tree(
                            original_tree,
                            ingroup_genes,
                            outgroup_genes,
                        )
                        outgroup_status = "ok"
                    else:
                        outgroup_status = "skip_no_valid_outgroup"
                    outgroup_report_rows.append(
                        (
                            tre_name,
                            parse_original_tree_id(tre_name),
                            str(len(ingroup_genes)),
                            str(len(outgroup_genes)),
                            ",".join(outgroup_genes),
                            outgroup_status,
                        )
                    )

                processed_lines.append(
                    tre_name + "\t" + serialize_tree_by_input_branch_length_style(
                        tree_to_write,
                        source_tree_path=tre_path,
                        fmt=0,
                    ) + "\n"
                )

                leaf_names = [leaf.name for leaf in Phylo_t_OG_L]
                tree_ortholog_trees.append(",".join(leaf_names))

            if tree_ortholog_trees:
                tsv_data[tre_ID] = tree_ortholog_trees

        sorted_lines = sorted(processed_lines, key=lambda x: len(x.split("\t")[1]), reverse=True)
        for line in sorted_lines:
            o.write(line)

    if synteny_blocks_path:
        with open("ortholog_synteny_report.tsv", "w", encoding="utf-8") as out:
            out.write(
                "tre_ID\ttre_name\tnum_genes\tsynteny_supported_genes\tbest_block\t"
                "best_block_gene_count\tsupport_ratio\tsynteny_filter\n"
            )
            for row in synteny_report_rows:
                out.write("\t".join(row) + "\n")

    if add_outgroup and outgroup_report_rows:
        with open("ortholog_outgroup_report.tsv", "w", encoding="utf-8") as out:
            out.write("tre_name\tgf_id\tnum_ingroup\tnum_outgroup\toutgroup_genes\tstatus\n")
            for row in outgroup_report_rows:
                out.write("\t".join(row) + "\n")

    if tsv_data:
        max_cols = max(len(trees) for trees in tsv_data.values())

        with open("ortholog_trees.tsv", "w", newline="", encoding="utf-8") as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')

            header = ["tre_ID"] + [f"OG_{i+1}" for i in range(max_cols)]
            writer.writerow(header)

            for tre_ID, trees in tsv_data.items():
                row = [tre_ID] + trees + [""] * (max_cols - len(trees))
                writer.writerow(row)


# ======================================================
# Section 5: CLI Entry Point
# ======================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Retrieve single-copy orthologs from gene family trees.",
    )
    parser.add_argument(
        "imap_file",
        help="Path to the imap file for gene ID transfer.",
    )
    parser.add_argument(
        "length_file",
        help="Path to the sequence length file.",
    )
    parser.add_argument(
        "gene_family_file",
        help="Path to the gene family list file.",
    )
    args = parser.parse_args()

    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic, _ = gene_id_transfer(args.imap_file)
    len_dic = read_and_return_dict(args.length_file)
    renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
    tre_dic = read_and_return_dict(args.gene_family_file)
    split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic, renamed_len_dic)
