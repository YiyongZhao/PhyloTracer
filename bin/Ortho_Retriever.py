"""
Ortholog retrieval and single-copy tree extraction for PhyloTracer analyses.

This module identifies offcut paralogous clades, separates ortholog groups,
extracts single-copy subtrees, and writes summary tables for downstream use.
"""

from __init__ import *

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
    for ev_seqs1 in paralogs_L:
        for ev_seqs2 in paralogs_L:
            if ev_seqs1.issubset(ev_seqs2):
                if ev_seqs1 not in paralogs_L:
                    paralogs_L.remove(ev_seqs1)
            elif ev_seqs2.issubset(ev_seqs1):
                if ev_seqs2 not in paralogs_L:
                    paralogs_L.remove(ev_seqs2)
    return paralogs_L


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


# ======================================================
# Section 4: Main Pipeline (Orchestration)
# ======================================================


def split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic, renamed_len_dic):
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

    o = open("ortho_retriever_summary.txt", "w")
    o.write("tre_name" + "\t" + "single_copy_tree" + "\n")

    processed_lines = []
    tsv_data = {}

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

        tree_ortholog_trees = []
        for clade in trees:
            tre_name = clade[0]
            Phylo_t_OG_L = clade[1]
            processed_lines.append(tre_name + "\t" + Phylo_t_OG_L.write() + "\n")

            leaf_names = [leaf.name for leaf in Phylo_t_OG_L]
            tree_ortholog_trees.append(",".join(leaf_names))

        if tree_ortholog_trees:
            tsv_data[tre_ID] = tree_ortholog_trees

    sorted_lines = sorted(processed_lines, key=lambda x: len(x.split("\t")[1]), reverse=True)
    for line in sorted_lines:
        o.write(line)
    o.close()

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
    gene2new_named_gene_dic, new_named_gene2gene_dic, voucher2taxa_dic = gene_id_transfer("imap")
    len_dic = read_and_return_dict("length")
    renamed_len_dic = rename_len_dic(len_dic, gene2new_named_gene_dic)
    tre_dic = read_and_return_dict("GF_list.txt")
    split_main(tre_dic, gene2new_named_gene_dic, new_named_gene2gene_dic, renamed_len_dic)
