"""
Ortholog retrieval and single-copy tree extraction for PhyloTracer analyses.

This module identifies offcut paralogous clades, separates ortholog groups,
extracts single-copy subtrees, and writes summary tables for downstream use.
"""

import re
import math
import hashlib
import os
import sqlite3
from collections import defaultdict
from pathlib import Path

from phylotracer import (
    find_tre_dup,
    is_rooted,
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


def _iter_synteny_pairs(synteny_blocks_path: str):
    """Yield explicit anchor pairs in file order, so the last duplicate wins."""
    with open(synteny_blocks_path, "r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Invalid synteny row at {synteny_blocks_path}:{line_number}")
            gene1 = parts[0]
            gene2 = parts[1]
            score = 1.0

            if len(parts) >= 3 and parts[1].lstrip("-").isdigit():
                gene2 = parts[2]
            elif len(parts) >= 3:
                try:
                    score = float(parts[2])
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid anchor score at {synteny_blocks_path}:{line_number}"
                    ) from exc
                if not math.isfinite(score) or score < 0:
                    raise ValueError(
                        f"Anchor score must be finite and nonnegative at "
                        f"{synteny_blocks_path}:{line_number}"
                    )

            if gene1 != gene2:
                yield gene1, gene2, score


def parse_synteny_blocks(synteny_blocks_path: str) -> dict[str, dict[str, float]]:
    """
    Parse explicit anchor pairs and their raw scores.

    Parameters
    ----------
    synteny_blocks_path : str
        Path to a raw block file where ``#`` starts a new block and each
        non-comment line stores a gene pair.

    Returns
    -------
    dict
        Mapping from each gene to partner genes and raw anchor scores. Pairs
        are stored in both directions; the last row wins for duplicate pairs.

    Assumptions
    -----------
    Each non-comment line contains either:
    1. a simplified two-column gene pair, or
    2. a WGDI-style row where the first and third columns are gene IDs.
    """
    gene_to_partners = defaultdict(dict)
    for gene1, gene2, score in _iter_synteny_pairs(synteny_blocks_path):
        gene_to_partners[gene1][gene2] = score
        gene_to_partners[gene2][gene1] = score

    return dict(gene_to_partners)


class SyntenySource:
    """Read one synteny file or index species-pair ``*.anchors`` files on demand."""

    def __init__(self, path: str, gene_to_species: dict[str, str]):
        source = Path(path)
        self.files_by_species_pair = {}
        self.covered_species_pairs = set()
        self._indexed = {}
        self._db = None
        self._single_scores = None

        code_to_species = {}
        ambiguous_codes = set()
        for gene, species in gene_to_species.items():
            if "_" in gene:
                code = gene.split("_", 1)[0]
                if code in ambiguous_codes:
                    continue
                if code in code_to_species and code_to_species[code] != species:
                    ambiguous_codes.add(code)
                    del code_to_species[code]
                else:
                    code_to_species[code] = species

        if source.is_dir():
            files = sorted(source.glob("*.anchors"))
            if not files:
                raise ValueError(f"No *.anchors files found in {source}")
            for file_path in files:
                tokens = file_path.stem.split(".")
                if len(tokens) != 2:
                    raise ValueError(f"Expected SpeciesA.SpeciesB.anchors: {file_path}")
                species_pair = tuple(sorted(code_to_species.get(token, token) for token in tokens))
                if species_pair in self.files_by_species_pair:
                    raise ValueError(f"Duplicate anchors files for species pair {species_pair}")
                self.files_by_species_pair[species_pair] = file_path
                self.covered_species_pairs.add(species_pair)
            cache_root = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
            cache_dir = cache_root / "phylotracer" / "synteny"
            cache_dir.mkdir(parents=True, exist_ok=True)
            directory_hash = hashlib.sha256(str(source.resolve()).encode("utf-8")).hexdigest()
            self._db = sqlite3.connect(cache_dir / f"{directory_hash}.sqlite3", timeout=60)
            self._db.execute("PRAGMA journal_mode=WAL")
            self._db.execute("PRAGMA synchronous=NORMAL")
            self._db.execute(
                "CREATE TABLE IF NOT EXISTS files ("
                "id INTEGER PRIMARY KEY, filename TEXT UNIQUE NOT NULL, "
                "size INTEGER NOT NULL, mtime_ns INTEGER NOT NULL, has_pairs INTEGER NOT NULL)"
            )
            self._db.execute(
                "CREATE TABLE IF NOT EXISTS pairs ("
                "file_id INTEGER NOT NULL, gene1 TEXT NOT NULL, gene2 TEXT NOT NULL, "
                "score REAL NOT NULL, PRIMARY KEY (file_id, gene1, gene2)) WITHOUT ROWID"
            )
        else:
            self._single_scores = parse_synteny_blocks(str(source))
            tokens = source.stem.split(".") if source.suffix == ".anchors" else []
            if len(tokens) == 2:
                self.covered_species_pairs.add(
                    tuple(sorted(code_to_species.get(token, token) for token in tokens))
                )
            else:
                for gene1, partners in self._single_scores.items():
                    species1 = gene_to_species.get(gene1)
                    if species1 is None:
                        continue
                    for gene2 in partners:
                        species2 = gene_to_species.get(gene2)
                        if species2 is not None and species1 != species2:
                            self.covered_species_pairs.add(tuple(sorted((species1, species2))))

    def _indexed_file(self, file_path: Path) -> tuple[int, bool]:
        """Index a file once per run and rebuild it if changed since a prior run."""
        filename = file_path.name
        if filename in self._indexed:
            return self._indexed[filename]
        stat = file_path.stat()
        row = self._db.execute(
            "SELECT id, size, mtime_ns, has_pairs FROM files WHERE filename = ?", (filename,)
        ).fetchone()
        if row and (row[1], row[2]) == (stat.st_size, stat.st_mtime_ns):
            result = (row[0], bool(row[3]))
            self._indexed[filename] = result
            return result

        with self._db:
            self._db.execute("BEGIN IMMEDIATE")
            row = self._db.execute(
                "SELECT id, size, mtime_ns, has_pairs FROM files WHERE filename = ?", (filename,)
            ).fetchone()
            if row and (row[1], row[2]) == (stat.st_size, stat.st_mtime_ns):
                result = (row[0], bool(row[3]))
                self._indexed[filename] = result
                return result
            if row:
                file_id = row[0]
                self._db.execute("DELETE FROM pairs WHERE file_id = ?", (file_id,))
            else:
                file_id = self._db.execute(
                    "INSERT INTO files (filename, size, mtime_ns, has_pairs) VALUES (?, ?, ?, 0)",
                    (filename, stat.st_size, stat.st_mtime_ns),
                ).lastrowid
            batch = []
            has_pairs = False
            for gene1, gene2, score in _iter_synteny_pairs(str(file_path)):
                has_pairs = True
                batch.append((file_id, min(gene1, gene2), max(gene1, gene2), score))
                if len(batch) == 1000:
                    self._db.executemany(
                        "INSERT OR REPLACE INTO pairs VALUES (?, ?, ?, ?)", batch
                    )
                    batch.clear()
            if batch:
                self._db.executemany("INSERT OR REPLACE INTO pairs VALUES (?, ?, ?, ?)", batch)
            current_stat = file_path.stat()
            if (current_stat.st_size, current_stat.st_mtime_ns) != (stat.st_size, stat.st_mtime_ns):
                raise RuntimeError(f"Synteny file changed while being indexed: {file_path}")
            self._db.execute(
                "UPDATE files SET size = ?, mtime_ns = ?, has_pairs = ? WHERE id = ?",
                (stat.st_size, stat.st_mtime_ns, int(has_pairs), file_id),
            )
        result = (file_id, has_pairs)
        self._indexed[filename] = result
        return result

    def close(self):
        if self._db is not None:
            self._db.close()
            self._db = None

    def pair_score(self, gene1: str, gene2: str, species1: str, species2: str) -> tuple[bool, float | None]:
        """Return whether the species pair was assayed and its optional anchor score."""
        species_pair = tuple(sorted((species1, species2)))
        if species_pair not in self.covered_species_pairs:
            return False, None

        if self._single_scores is not None:
            scores = self._single_scores
            if not scores:
                return False, None
            return True, scores.get(gene1, {}).get(gene2)

        file_id, has_pairs = self._indexed_file(self.files_by_species_pair[species_pair])
        if not has_pairs:
            return False, None
        row = self._db.execute(
            "SELECT score FROM pairs WHERE file_id = ? AND gene1 = ? AND gene2 = ?",
            (file_id, min(gene1, gene2), max(gene1, gene2)),
        ).fetchone()
        return True, row[0] if row else None


def filter_trees_by_synteny(
    trees: list,
    gene_to_partners: SyntenySource | dict,
    gene_to_species: dict[str, str],
    synteny_threshold: float = 0.2,
) -> tuple[list, list]:
    """
    Filter candidate ortholog trees by raw synteny block support.

    Parameters
    ----------
    trees : list
        List of ``(tree_name, tree_object)`` tuples.
    gene_to_partners : SyntenySource or dict
        Synteny source or a direct gene-to-partner mapping.
    gene_to_species : dict
        Species mapping for every candidate gene. Each candidate must contain
        at most one gene per species.
    synteny_threshold : float
        Minimum matched fraction of all candidate cross-species pairs. A pair
        missing from the synteny data counts as unmatched.

    Returns
    -------
    tuple[list, list]
        Filtered trees and per-tree synteny report rows.

    Assumptions
    -----------
    A valid single-copy tree with n tips has n * (n - 1) / 2 pairs. The rate is
    ``matched_pairs / total_pairs``. Missing genes, pairs, or species-pair
    files never reduce the denominator. Raw anchor scores are reported only as
    a separate mean over matched pairs.
    """
    if not 0.0 <= synteny_threshold <= 1.0:
        raise ValueError("synteny_threshold must be between 0 and 1")

    filtered_trees = []
    report_rows = []

    for tree_name, phylo_tree in trees:
        leaf_names = [leaf.name for leaf in phylo_tree]
        if len(leaf_names) != len(set(leaf_names)):
            raise ValueError(f"Duplicate gene IDs in candidate tree {tree_name}")
        missing_species = [gene for gene in leaf_names if not gene_to_species.get(gene)]
        if missing_species:
            raise ValueError(
                f"Missing species mapping for candidate tree {tree_name}: "
                + ", ".join(missing_species)
            )
        species_names = [gene_to_species[gene] for gene in leaf_names]
        if len(species_names) != len(set(species_names)):
            raise ValueError(f"Candidate tree {tree_name} is not single-copy: duplicate species")
        total_pairs = len(leaf_names) * (len(leaf_names) - 1) // 2
        matched_pairs = 0
        pairs_without_synteny_data = 0
        total_score = 0.0
        for index, gene1 in enumerate(leaf_names):
            species1 = gene_to_species[gene1]
            for gene2 in leaf_names[index + 1:]:
                species2 = gene_to_species[gene2]
                if isinstance(gene_to_partners, SyntenySource):
                    covered, score = gene_to_partners.pair_score(
                        gene1, gene2, species1, species2
                    )
                    if not covered:
                        pairs_without_synteny_data += 1
                else:
                    partners = gene_to_partners.get(gene1, {})
                    reverse_partners = gene_to_partners.get(gene2, {})
                    score = (
                        partners.get(gene2) if isinstance(partners, dict)
                        else 1.0 if gene2 in partners else None
                    )
                    if score is None:
                        score = (
                            reverse_partners.get(gene1) if isinstance(reverse_partners, dict)
                            else 1.0 if gene1 in reverse_partners else None
                        )
                if score is not None:
                    matched_pairs += 1
                    total_score += score
        unmatched_pairs = total_pairs - matched_pairs
        synteny_rate = matched_pairs / total_pairs if total_pairs else 0.0
        avg_matched_score = total_score / matched_pairs if matched_pairs else 0.0
        passes = total_pairs > 0 and synteny_rate >= synteny_threshold
        if passes:
            filtered_trees.append((tree_name, phylo_tree))

        report_rows.append(
            (
                tree_name,
                str(len(leaf_names)),
                str(total_pairs),
                str(matched_pairs),
                str(unmatched_pairs),
                str(pairs_without_synteny_data),
                f"{synteny_rate:.6f}",
                f"{avg_matched_score:.2f}",
                f"{synteny_threshold:.4f}",
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
    synteny_threshold: float = 0.2,
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
    synteny_threshold : float
        Minimum cross-species gene-pair rate required to retain a subtree when
        a synteny file is supplied.

    Returns
    -------
    None

    Assumptions
    -----------
    Input trees can be rooted and pruned without identifier conflicts.
    """
    import csv

    if synteny_blocks_path and gene_to_species is None:
        raise ValueError("gene_to_species is required when synteny filtering is enabled")
    if synteny_blocks_path and not 0.0 <= synteny_threshold <= 1.0:
        raise ValueError("synteny_threshold must be between 0 and 1")

    gene_to_partners = None
    if synteny_blocks_path:
        gene_to_partners = SyntenySource(synteny_blocks_path, gene_to_species)

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
                trees, tree_report_rows = filter_trees_by_synteny(
                    trees,
                    gene_to_partners,
                    gene_to_species=gene_to_species,
                    synteny_threshold=synteny_threshold,
                )
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
                "tre_ID\ttre_name\tnum_genes\ttotal_pairs\tmatched_pairs\t"
                "unmatched_pairs\tpairs_without_synteny_data\tsynteny_rate\t"
                "mean_matched_anchor_score\tthreshold\tsynteny_filter\n"
            )
            for row in synteny_report_rows:
                out.write("\t".join(row) + "\n")

    if add_outgroup:
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
    elif synteny_blocks_path:
        with open("ortholog_trees.tsv", "w", newline="", encoding="utf-8") as tsvfile:
            csv.writer(tsvfile, delimiter='\t').writerow(["tre_ID"])

    if isinstance(gene_to_partners, SyntenySource):
        gene_to_partners.close()


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
