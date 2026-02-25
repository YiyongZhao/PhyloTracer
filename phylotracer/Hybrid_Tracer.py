"""
Hybridization tracing and HyDe integration for the PhyloTracer pipeline.

This module identifies candidate duplication clades, builds gene matrices,
runs HyDe analyses, and writes filtered hybridization results.
"""

import logging
import os
import subprocess
import time

logger = logging.getLogger(__name__)
from collections import Counter, defaultdict
from signal import SIGUSR2

import numpy as np
import pandas as pd
try:
    import phyde as hd
except ImportError:
    hd = None
from Bio import SeqIO
from ete3 import PhyloTree
from tqdm import tqdm

from phylotracer import (
    annotate_gene_tree,
    find_dup_node,
    gene_id_transfer,
    get_species_list,
    get_species_set,
    num_sptree,
    num_tre_node,
    read_and_return_dict,
    read_phylo_tree,
    rename_input_tre,
)
from phylotracer.GD_Detector import get_model as detector_get_model, normalize_model

# ======================================================
# Section 1: Species Tree and Sequence Utilities
# ======================================================


def process_start_node(file_path: str, sptree: object) -> str:
    """
    Determine the common ancestor of species listed in a file.

    Parameters
    ----------
    file_path : str
        Path to a file containing species names, one per line.
    sptree : object
        Species tree object used to compute the common ancestor.

    Returns
    -------
    str
        Name of the common ancestor node, or None if unavailable.

    Assumptions
    -----------
    Species names in the file match those in the species tree.
    """
    species_list = []
    try:
        with open(file_path, "r") as f:
            for line in f:
                species_list.append(line.strip())
    except FileNotFoundError:
        logger.error("Species list file not found at %s", file_path)
        return None
    except Exception as e:
        logger.error("Error reading species list file: %s", e)
        return None

    if not species_list:
        logger.warning("Species list file is empty.")
        return None

    try:
        common_ancestor = sptree.get_common_ancestor(species_list)
        return common_ancestor.name
    except Exception as e:
        logger.error("Error finding common ancestor in species tree: %s", e)
        return None


def create_fasta_dict(fasta_file, gene2new_named_gene_dic):
    """
    Parse a FASTA file into a dictionary keyed by renamed gene IDs.

    Parameters
    ----------
    fasta_file : object
        Path to a FASTA file.
    gene2new_named_gene_dic : dict
        Mapping from original gene identifiers to renamed identifiers.

    Returns
    -------
    dict
        Mapping from renamed gene identifiers to sequence strings.

    Assumptions
    -----------
    All record IDs are present in ``gene2new_named_gene_dic``.
    """
    fasta_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        fasta_dict[gene2new_named_gene_dic[record.id]] = str(record.seq)
    return fasta_dict


def get_outsps_sort_lst(map_t, sptree):
    """
    Sort outgroup species by topological distance to a mapped node.

    Parameters
    ----------
    map_t : object
        Species tree node defining the mapped clade.
    sptree : object
        Species tree object to traverse.

    Returns
    -------
    list
        List of species names sorted by distance to the mapped node.

    Assumptions
    -----------
    Distances are computed using topology-only mode.
    """
    sps_depth = {}
    for gene in sptree:
        if gene.name in get_species_list(map_t):
            continue
        sps_depth[gene.name] = gene.get_distance(map_t, topology_only=True)
    sorted_out = sorted(sps_depth.items(), key=lambda item: item[1])
    return [sps for sps, _ in sorted_out]


def get_dynamic_basal_set_for_gd(gene_tree_species_set: set, species_tree_root: object) -> set:
    """
    Compute a dynamic outgroup-candidate species set using Phylo_Rooter logic.

    Parameters
    ----------
    gene_tree_species_set : set
        Species covered by the current gene tree.
    species_tree_root : object
        Root of the species tree.

    Returns
    -------
    set
        Candidate outgroup species set for the current gene tree.

    Assumptions
    -----------
    Species tree is rooted and generally bifurcating.
    """
    if species_tree_root.is_leaf():
        return {species_tree_root.name}

    children = species_tree_root.get_children()
    clade_sets = [set(child.get_leaf_names()) for child in children]
    present_flags = [not gene_tree_species_set.isdisjoint(c_set) for c_set in clade_sets]

    if sum(present_flags) >= 2:
        valid_children = [
            (children[i], len(clade_sets[i])) for i, has_gene in enumerate(present_flags) if has_gene
        ]
        valid_children.sort(key=lambda x: x[1])
        return set(valid_children[0][0].get_leaf_names())

    if sum(present_flags) == 1:
        for i, has_gene in enumerate(present_flags):
            if has_gene:
                return get_dynamic_basal_set_for_gd(gene_tree_species_set, children[i])

    return set()


def get_outsps_min_distance_gene(node, sps):
    """
    Identify the closest gene in a node matching a species prefix.

    Parameters
    ----------
    node : object
        Tree node containing candidate genes.
    sps : object
        Species prefix to match.

    Returns
    -------
    object
        Gene name with minimal distance, or None if absent.

    Assumptions
    -----------
    Gene names include species identifiers as prefixes.
    """
    dic = {i.name: node.get_distance(i) for i in node if sps in i.name}
    if not dic:
        return None
    return min(dic, key=dic.get)


def get_sister_nodes(node):
    """
    Collect sister nodes along the path to the root.

    Parameters
    ----------
    node : object
        Tree node from which to traverse upward.

    Returns
    -------
    list
        List of sister nodes encountered on the path to the root.

    Assumptions
    -----------
    Tree is connected and node has a valid ancestor chain.
    """
    sister_nodes = []
    current_node = node
    while current_node and not current_node.is_root():
        sisters = current_node.get_sisters()
        if sisters:
            sister_nodes.append(sisters[0])
        current_node = current_node.up
    return sister_nodes


def get_target_gene(sisters, sps):
    """
    Find a target gene within sister clades for a given species.

    Parameters
    ----------
    sisters : list
        List of sister nodes to search.
    sps : object
        Species prefix to match.

    Returns
    -------
    object
        Closest matching gene name, or None if not found.

    Assumptions
    -----------
    Sister nodes are valid ETE nodes with leaf names.
    """
    for sis_node in sisters:
        if sps in get_species_set(sis_node):
            return get_outsps_min_distance_gene(sis_node, sps)
    return None


def get_outgroup_species_for_gd_node(gd_node_name, sptree):
    """
    Select the nearest sister-branch species as fixed outgroup for a GD node.

    Parameters
    ----------
    gd_node_name : str
        Species-tree node name to which the GD group is mapped.
    sptree : object
        Species tree used for reconciliation.

    Returns
    -------
    str or None
        Selected outgroup species name, or None if unavailable.
    """
    map_t = sptree & gd_node_name
    sister_nodes = map_t.get_sisters()
    if not sister_nodes:
        return None

    sister_species = sister_nodes[0].get_leaf_names()
    if not sister_species:
        return None

    species_distances = []
    for species in sister_species:
        species_node = sptree & species
        distance = map_t.get_distance(species_node, topology_only=True)
        species_distances.append((distance, species))

    species_distances.sort(key=lambda x: (x[0], x[1]))
    return species_distances[0][1]


def get_outgroup_gene(gd_clade, outgroup_species):
    """
    Identify the closest outgroup gene of a fixed outgroup species.

    Parameters
    ----------
    gd_clade : object
        Duplication clade in the gene tree.
    outgroup_species : str
        Fixed outgroup species for the current GD node.

    Returns
    -------
    str or None
        Closest outgroup gene name, or None if absent.

    Assumptions
    -----------
    Candidate outgroup genes are searched only in sister lineages of the GD
    clade along the path to the root.
    """
    sisters = get_sister_nodes(gd_clade)
    best_gene = None
    best_distance = None

    for sister_node in sisters:
        for leaf in sister_node.get_leaves():
            species = leaf.name.split("_")[0]
            if species != outgroup_species:
                continue
            distance = gd_clade.get_distance(leaf)
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_gene = leaf.name

    return best_gene


# ======================================================
# Section 2: Duplication Model Assignment
# ======================================================


def get_model(clade, sptree):
    """
    Assign GD type using the project-standard detector definitions.

    Parameters
    ----------
    clade : object
        Gene clade with mapped species-tree label.
    sptree : object
        Species tree used for reconciliation.

    Returns
    -------
    str
        Canonical GD type (AABB/AXBB/AABX/Complex).
    """
    return normalize_model(detector_get_model(clade, sptree))


# ======================================================
# Section 3: GD Counting and Filtering
# ======================================================


def count_elements_in_lists(data):
    """
    Count GD model labels across clades using canonical detector types.

    Parameters
    ----------
    data : dict
        Mapping from node to list of model labels.

    Returns
    -------
    dict
        Mapping from node to Counter of canonical model types.

    Assumptions
    -----------
    Input values are canonical GD types from ``GD_Detector.normalize_model``.
    """
    return {key: Counter(value) for key, value in data.items()}


def get_gd_count_dic_and_gd_type_dic(
    tre_dic,
    gene2new_named_gene_dic,
    rename_sptree,
):
    """
    Collect duplication counts and model types across gene trees.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to file paths.
    gene2new_named_gene_dic : dict
        Mapping from gene identifiers to renamed identifiers.
    rename_sptree : object
        Renamed species tree used for mapping.
    Returns
    -------
    tuple
        (gd_count_dic, gd_type_dic) dictionaries.

    Assumptions
    -----------
    Duplication nodes are detected by ``find_dup_node`` with default support.
    """

    gd_num = 0
    gd_count_dic = {}
    gd_type_dic = {}
    for tre_id, tre_path in tre_dic.items():
        t = read_phylo_tree(tre_path)
        t1 = rename_input_tre(t, gene2new_named_gene_dic)
        annotate_gene_tree(t1, rename_sptree)
        gds = find_dup_node(t1, rename_sptree,50,0,2,0,1)

        for gd in gds:
            mapped_node = rename_sptree & gd.map
            if len(get_species_set(mapped_node)) < 2:
                continue
            if len(mapped_node.get_children()) < 2:
                continue
            try:
                type_str = get_model(gd, rename_sptree)
            except Exception:
                continue
            if gd.map in gd_count_dic:
                gd_count_dic[gd.map].append((tre_id + "-" + str(gd_num), gd))
                gd_type_dic[gd.map].append(type_str)
            else:
                gd_count_dic[gd.map] = [(tre_id + "-" + str(gd_num), gd)]
                gd_type_dic[gd.map] = [type_str]
            gd_num += 1
    return gd_count_dic, gd_type_dic


def get_process_gd_clade(gd_type_dic, gd_count_dic):
    """
    Filter GD clades based on model composition thresholds.

    Parameters
    ----------
    gd_type_dic : dict
        Mapping from clade to model type counters.
    gd_count_dic : dict
        Mapping from clade to GD event lists.

    Returns
    -------
    list
        List of filtered GD clades with event lists.

    Assumptions
    -----------
    Thresholds follow the original pipeline heuristics.
    """
    gd_clade = []
    for k, v in gd_type_dic.items():
        asymmetric_num = v["AXBB"] + v["AABX"]
        total_num = v["AXBB"] + v["AABX"] + v["AABB"]
        if total_num != 0:
            type_ratio = asymmetric_num / total_num
            if total_num >= 10 and type_ratio >= 0.1:
                gd_clade.append((k, gd_count_dic[k]))
    return gd_clade


# ======================================================
# Section 4: HyDe Output Formatting
# ======================================================


def write_out(out, triple, outfile):
    """
    Write a HyDe result record to an output file.

    Parameters
    ----------
    out : dict
        HyDe result dictionary.
    triple : tuple
        Taxon triple tuple.
    outfile : object
        Open file handle for output.

    Returns
    -------
    None

    Assumptions
    -----------
    HyDe output dictionary contains required keys.
    """
    print(triple[0], "\t", triple[1], "\t", triple[2], "\t", sep="", end="", file=outfile)
    print(out["Zscore"], "\t", sep="", end="", file=outfile)
    print(out["Pvalue"], "\t", sep="", end="", file=outfile)
    print(out["Gamma"], "\t", sep="", end="", file=outfile)
    print(out["AAAA"], "\t", sep="", end="", file=outfile)
    print(out["AAAB"], "\t", sep="", end="", file=outfile)
    print(out["AABA"], "\t", sep="", end="", file=outfile)
    print(out["AABB"], "\t", sep="", end="", file=outfile)
    print(out["AABC"], "\t", sep="", end="", file=outfile)
    print(out["ABAA"], "\t", sep="", end="", file=outfile)
    print(out["ABAB"], "\t", sep="", end="", file=outfile)
    print(out["ABAC"], "\t", sep="", end="", file=outfile)
    print(out["ABBA"], "\t", sep="", end="", file=outfile)
    print(out["BAAA"], "\t", sep="", end="", file=outfile)
    print(out["ABBC"], "\t", sep="", end="", file=outfile)
    print(out["CABC"], "\t", sep="", end="", file=outfile)
    print(out["BACA"], "\t", sep="", end="", file=outfile)
    print(out["BCAA"], "\t", sep="", end="", file=outfile)
    print(out["ABCD"], "\n", sep="", end="", file=outfile)


# ======================================================
# Section 5: HyDe Processing Pipeline
# ======================================================


def hyde_main(
    tre_dic,
    seq_path_dic,
    rename_sptree,
    gene2new_named_gene_dic,
    voucher2taxa_dic,
    target_node=None,
    gd_group: int = 1,
):
    """
    Run HyDe analyses on filtered duplication clades.

    Parameters
    ----------
    tre_dic : dict
        Mapping from tree IDs to file paths.
    seq_path_dic : dict
        Mapping from tree IDs to FASTA file paths.
    rename_sptree : object
        Renamed species tree used for mapping.
    gene2new_named_gene_dic : dict
        Mapping from gene identifiers to renamed identifiers.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.
    target_node : object, optional
        Specific species-tree node to focus on.
    gd_group : int, optional
        Number of parallel processing groups.
    Returns
    -------
    None

    Assumptions
    -----------
    HyDe results are written to ``hyde_out.txt`` and ``hyde_filtered_out.txt``.
    """
    hyde_tuple_lst = []
    gd_count_dic, gd_type_dic = get_gd_count_dic_and_gd_type_dic(
        tre_dic,
        gene2new_named_gene_dic,
        rename_sptree,
    )
    data = count_elements_in_lists(gd_type_dic)
    gd_clades = get_process_gd_clade(data, gd_count_dic)

    for gd in gd_clades:
        gd_name, gds = gd
        if target_node and gd_name != target_node:
            continue

        logger.info("%s is processing", gd_name)
        logger.info("Total gene duplication events in dataset: %d", len(gds))
        logger.info("Number of parallel processing groups: %d", gd_group)

        split_gds = np.array_split(gds, gd_group)
        fixed_outgroup_species = get_outgroup_species_for_gd_node(gd_name, rename_sptree)
        if fixed_outgroup_species is None:
            logger.warning(
                "[Hybrid_Tracer][Skip] GD node %s has no sister-branch "
                "outgroup candidates. All events are skipped.", gd_name
            )
            continue
        fixed_outgroup_original = voucher2taxa_dic.get(
            fixed_outgroup_species,
            fixed_outgroup_species,
        )
        logger.info(
            "[Hybrid_Tracer] GD node %s uses fixed outgroup species: %s",
            gd_name, fixed_outgroup_original
        )

        logger.info("Gene duplication events partitioned into %d processing batches", len(split_gds))
        for i, sub_group in enumerate(split_gds):
            logger.info("Batch %d: %d gene duplication events", i+1, len(sub_group))
        for sub_group in split_gds:
            hyde_result_lst = process_gd_group(
                sub_group,
                rename_sptree,
                gd_name,
                seq_path_dic,
                gene2new_named_gene_dic,
                voucher2taxa_dic,
                target_node,
                fixed_outgroup_species,
            )
            hyde_tuple_lst.extend(hyde_result_lst)
    header = (
        "P1\tHybrid\tP2\tZscore\tPvalue\tGamma\tAAAA\tAAAB\tAABA\tAABB\tAABC\tABAA\t"
        "ABAB\tABAC\tABBA\tBAAA\tABBC\tCABC\tBACA\tBCAA\tABCD\n"
    )
    with open("hyde_out.txt", "w") as out_file, open(
        "hyde_filtered_out.txt", "w"
    ) as filtered_out_file:
        print(header, end="", file=out_file)
        print(header, end="", file=filtered_out_file)

        for t, res, t_num in hyde_tuple_lst:
            write_out(res, t, out_file)
            if is_filtered(res, t_num):
                write_out(res, t, filtered_out_file)


def process_gd_group(
    gds,
    rename_sptree,
    gd_name,
    seq_path_dic,
    gene2new_named_gene_dic,
    voucher2taxa_dic,
    target_node,
    fixed_outgroup_species=None,
):
    """
    Process a subset of GD events for HyDe analysis.

    Parameters
    ----------
    gds : list
        List of GD events for processing.
    rename_sptree : object
        Renamed species tree.
    gd_name : object
        GD clade name.
    seq_path_dic : dict
        Mapping from tree IDs to FASTA file paths.
    gene2new_named_gene_dic : dict
        Mapping from gene identifiers to renamed identifiers.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.
    target_node : object
        Target species-tree node used for filtering.
    fixed_outgroup_species : str, optional
        Fixed outgroup species used in strict mode.
    Returns
    -------
    list
        List of HyDe results for the processed group.

    Assumptions
    -----------
    Outgroup selection uses ``get_outgroup_gene`` heuristics.
    """
    target_node_name = target_node if target_node else gd_name
    target_clade = rename_sptree & target_node_name
    all_matrices = []
    col_counter = 1
    seq_dic_all = {}
    skipped_no_outgroup = 0

    for gd_clade_set in gds:
        gdid = gd_clade_set[0]
        gd_clade = gd_clade_set[1]

        outfile = gdid
        tre_id1 = outfile.split("-")[0]

        try:
            seq_dic = create_fasta_dict(seq_path_dic[tre_id1], gene2new_named_gene_dic)
            seq_dic_all.update(seq_dic)
        except KeyError as e:
            logger.warning("KeyError encountered for GD %s: %s. Skipping this GD.", gdid, e)
            continue

        outgroup_species = fixed_outgroup_species
        outgroup_gene = get_outgroup_gene(gd_clade, outgroup_species)
        if outgroup_gene is None:
            skipped_no_outgroup += 1
            gd_map = rename_sptree.get_common_ancestor(get_species_set(gd_clade))
            if outgroup_species is None:
                expected_outgroup_species = "None"
            else:
                expected_outgroup_species = voucher2taxa_dic.get(
                    outgroup_species,
                    outgroup_species,
                )
            logger.warning(
                "[Hybrid_Tracer][Skip] GD event %s (map=%s) has no valid outgroup gene "
                "from expected outgroup species %s.", gdid, gd_map.name, expected_outgroup_species
            )
            continue
        one_gd_matrix = process_one_gd(gd_clade, outgroup_gene)
        num_cols = one_gd_matrix.shape[1]
        new_col_names = [f"col{col_counter + i}" for i in range(num_cols)]
        one_gd_matrix.columns = new_col_names

        col_counter += num_cols

        all_matrices.append(one_gd_matrix)

    if skipped_no_outgroup > 0:
        logger.info(
            "[Hybrid_Tracer][Summary] GD node %s: skipped %d events "
            "due to missing valid outgroup genes.", gd_name, skipped_no_outgroup
        )

    if not all_matrices:
        logger.warning("No valid matrices found for GD group %s. Returning empty result.", gd_name)
        return []

    merged_matrix = pd.concat(all_matrices, axis=1)
    merged_matrix = merged_matrix.fillna("-")

    clean_matrix = clean_matrix_by_dash_count(merged_matrix)
    hyde_result_lst = run_hyde_from_matrix_integrated(
        clean_matrix,
        seq_dic_all,
        voucher2taxa_dic,
        target_clade,
    )
    return hyde_result_lst


def run_hyde_from_matrix_integrated(
    clean_matrix: pd.DataFrame,
    seq_dic: dict,
    voucher2taxa_dic: dict,
    clade: object,
    trim: bool = True,
):
    """
    Run HyDe analysis from a cleaned gene matrix.

    Parameters
    ----------
    clean_matrix : pd.DataFrame
        Cleaned gene matrix with species as rows.
    seq_dic : dict
        Mapping from gene identifiers to sequences.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.
    clade : object
        Target clade defining in-group species.
    trim : bool, optional
        Whether to trim columns of all gaps.

    Returns
    -------
    list
        List of HyDe result tuples.

    Assumptions
    -----------
    Temporary files ``temp.phy`` and ``temp.imap`` can be created.
    """
    if clean_matrix is None or clean_matrix.empty or clean_matrix.shape[1] == 0:
        return []
    if not seq_dic:
        return []

    target_sps = clade.get_leaf_names()

    import tempfile

    tmp_phy = tempfile.NamedTemporaryFile(mode='w', suffix='.phy', delete=False)
    tmp_imap = tempfile.NamedTemporaryFile(mode='w', suffix='.imap', delete=False)
    tmp_phy_path = tmp_phy.name
    tmp_imap_path = tmp_imap.name
    tmp_phy.close()
    tmp_imap.close()

    matrix_to_phy(clean_matrix, seq_dic, voucher2taxa_dic, tmp_phy_path, trim=trim)

    imap = {}
    with open(tmp_imap_path, "w") as imap_file:
        for species in clean_matrix.index:
            sp_0 = voucher2taxa_dic[species]
            if species not in target_sps:
                imap[sp_0] = "out"
                imap_file.write(f"{sp_0}\tout\n")
            else:
                imap[sp_0] = sp_0
                imap_file.write(f"{sp_0}\t{sp_0}\n")

    with open(tmp_phy_path, "r") as f:
        first_line = f.readline().strip()
        num_species, max_length = map(int, first_line.split())

    hyde_result_lst = []
    sps_num = len(imap.keys())
    taxa_num = len(set(imap.values()))

    try:
        dat = hd.HydeData(tmp_phy_path, tmp_imap_path, "out", sps_num, taxa_num, max_length, quiet=True)
        res = dat.list_triples()

        for t in res:
            p1, hyb, p2 = t
            result = dat.test_triple(p1, hyb, p2)
            combined_element = (t, result, len(res))
            hyde_result_lst.append(combined_element)
    finally:
        os.remove(tmp_phy_path)
        os.remove(tmp_imap_path)

    return hyde_result_lst


def clean_matrix_by_dash_count(matrix: pd.DataFrame, max_allowed_dashes: int = 1):
    """
    Remove duplicate columns and columns with excessive missing data.

    Parameters
    ----------
    matrix : pd.DataFrame
        Input gene matrix.
    max_allowed_dashes : int, optional
        Maximum allowed dashes per column.

    Returns
    -------
    pd.DataFrame
        Cleaned gene matrix.

    Assumptions
    -----------
    Dash characters represent missing data.
    """
    matrix = matrix.T.drop_duplicates().T
    dash_counts = (matrix == "-").sum(axis=0)
    matrix = matrix.loc[:, dash_counts <= max_allowed_dashes]

    matrix.columns = [f"col{i+1}" for i in range(matrix.shape[1])]

    return matrix


def matrix_to_phy(
    clean_matrix: pd.DataFrame,
    seq_dic: dict,
    voucher2taxa_dic,
    output_file: str = "output.phy",
    trim: bool = True,
):
    """
    Convert a gene matrix to PHYLIP format.

    Parameters
    ----------
    clean_matrix : pd.DataFrame
        Cleaned gene matrix with species as rows.
    seq_dic : dict
        Mapping from gene identifiers to sequences.
    voucher2taxa_dic : dict
        Mapping from voucher identifiers to taxa labels.
    output_file : str, optional
        Output PHYLIP file path.
    trim : bool, optional
        Whether to trim columns of all gaps.

    Returns
    -------
    None

    Assumptions
    -----------
    Gene identifiers in the matrix correspond to keys in ``seq_dic``.
    """
    species_seqs = defaultdict(str)

    for col in clean_matrix.columns:
        for species, gene_id in clean_matrix[col].items():
            sp = voucher2taxa_dic[species]
            if gene_id != "-" and gene_id in seq_dic:
                seq = seq_dic[gene_id]
            else:
                seq = "-" * max(len(s) for s in seq_dic.values())
            species_seqs[sp] += seq

    if trim:
        def trimal_matrix(concat_matrix: dict):
            """
            Remove columns consisting entirely of gaps.

            Parameters
            ----------
            concat_matrix : dict
                Mapping from species name to concatenated sequence.

            Returns
            -------
            tuple
                (filtered_matrix, new_length).

            Assumptions
            -----------
            Sequences can be padded to a common length with gaps.
            """
            max_len = max(len(seq) for seq in concat_matrix.values())

            for sp in concat_matrix:
                seq = concat_matrix[sp]
                if len(seq) < max_len:
                    concat_matrix[sp] = seq + "-" * (max_len - len(seq))

            valid_cols = [
                i
                for i in range(max_len)
                if any(concat_matrix[sp][i] != "-" for sp in concat_matrix)
            ]

            filtered_matrix = {
                sp: "".join(concat_matrix[sp][i] for i in valid_cols)
                for sp in concat_matrix
            }

            return filtered_matrix, len(valid_cols)

        species_seqs, trimmed_len = trimal_matrix(species_seqs)
    else:
        trimmed_len = len(next(iter(species_seqs.values())))

    num_species = len(species_seqs)
    with open(output_file, "w") as f:
        f.write(f"{num_species} {trimmed_len}\n")
        for species, seq in species_seqs.items():
            name = species[:10].ljust(10)
            f.write(f"{name} {seq}\n")


def process_one_gd(gd_clade: object, outgroup_gene: str):
    """
    Build a gene matrix for a single GD clade and outgroup gene.

    Parameters
    ----------
    gd_clade : object
        Duplication clade to process.
    outgroup_gene : str
        Outgroup gene identifier.

    Returns
    -------
    pd.DataFrame
        Gene matrix for the GD clade.

    Assumptions
    -----------
    Gene identifiers encode species prefixes separated by underscores.
    """
    result = []

    def is_dup_node(node):
        genes = get_species_list(node)
        species = get_species_set(node)
        return len(genes) != len(species)

    def traverse(node):
        if is_dup_node(node):
            for child in getattr(node, "children", []):
                traverse(child)
        else:
            genes = node.get_leaf_names() if hasattr(node, "get_leaf_names") else []
            tuple_genes = tuple(sorted(genes + [outgroup_gene], key=lambda x: x.split("_")[0]))
            result.append(tuple_genes)

    traverse(gd_clade)

    if not result:
        return pd.DataFrame()

    all_species = set()
    for triple in result:
        for gene in triple:
            species = gene.split("_")[0]
            all_species.add(species)

    all_species = sorted(list(all_species))

    matrix_data = []
    for i, triple in enumerate(result):
        column_data = []
        for species in all_species:
            genes_in_species = [gene for gene in triple if gene.split("_")[0] == species]
            if genes_in_species:
                column_data.append(";".join(genes_in_species))
            else:
                column_data.append("-")
        matrix_data.append(column_data)

    matrix_df = pd.DataFrame(matrix_data).T
    matrix_df.index = all_species
    matrix_df.columns = [f"col_{i+1}" for i in range(len(result))]

    return matrix_df


def is_filtered(res, t_num):
    """
    Check whether a HyDe result passes filtering thresholds.

    Parameters
    ----------
    res : dict
        HyDe result dictionary.
    t_num : int
        Number of tested triples.

    Returns
    -------
    bool
        True if the result passes filtering criteria.

    Assumptions
    -----------
    Bonferroni correction is applied to the p-value threshold.
    """
    return (
        res["Pvalue"] < (0.05 / t_num)
        and abs(res["Zscore"]) != 99999.9
        and 0.0 < res["Gamma"] < 1.0
    )
