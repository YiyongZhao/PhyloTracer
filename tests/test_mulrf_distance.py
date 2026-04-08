from __future__ import annotations

from ete3 import Tree

from phylotracer.MulRF_Distance import (
    compute_mulrf,
    compute_mulrf_between_gene_trees,
    compute_mulrf_gene_vs_species,
)


def _mk_tree(newick: str) -> Tree:
    return Tree(newick)


def _old_species_set_splits(tree: Tree, sp_map: dict[str, str]) -> set[frozenset[str]]:
    """
    Legacy-style approximation used only for test contrast:
    - Collapse descendants to plain species sets (ignore copy counts)
    - Store unique split sets only (ignore multiplicity)
    """
    total_species = {sp_map[leaf.name] for leaf in tree.get_leaves()}
    out: set[frozenset[str]] = set()
    for node in tree.traverse("postorder"):
        if node.is_leaf() or node.is_root():
            continue
        side = {sp_map[leaf.name] for leaf in node.get_leaves()}
        other = total_species - side
        if not side or not other:
            continue
        side_key = tuple(sorted(side))
        other_key = tuple(sorted(other))
        canonical = frozenset(side if side_key <= other_key else other)
        out.add(canonical)
    return out


def _old_species_collapse_distance(
    tree1: Tree,
    tree2: Tree,
    map1: dict[str, str],
    map2: dict[str, str],
) -> int:
    s1 = _old_species_set_splits(tree1, map1)
    s2 = _old_species_set_splits(tree2, map2)
    return len(s1.symmetric_difference(s2))


def test_single_copy_identical_gene_vs_species_zero_distance():
    # Identical single-copy topology should have zero conflict.
    gene_tree = _mk_tree("((A_g1,B_g1),(C_g1,D_g1));")
    species_tree = _mk_tree("((A,B),(C,D));")
    gmap = {"A_g1": "A", "B_g1": "B", "C_g1": "C", "D_g1": "D"}

    res = compute_mulrf(gene_tree, species_tree, gene2sp_map=gmap)
    assert res["mulrf"] == 0
    assert res["normalized_mulrf"] == 0.0


def test_single_copy_conflicting_gene_vs_species_nonzero_distance():
    # Conflicting single-copy topology should have positive distance.
    gene_same = _mk_tree("((A_g1,B_g1),(C_g1,D_g1));")
    gene_diff = _mk_tree("((A_g1,C_g1),(B_g1,D_g1));")
    species_tree = _mk_tree("((A,B),(C,D));")
    gmap = {"A_g1": "A", "B_g1": "B", "C_g1": "C", "D_g1": "D"}

    res_same = compute_mulrf(gene_same, species_tree, gene2sp_map=gmap)
    res_diff = compute_mulrf(gene_diff, species_tree, gene2sp_map=gmap)

    assert res_diff["mulrf"] > 0
    assert res_diff["shared_bipartitions"] < res_same["shared_bipartitions"]


def test_multi_copy_gene_vs_species_computable_and_counts_correct():
    # Multi-copy gene tree vs species tree should be computable with correct shared species.
    gene_tree = _mk_tree("(((A_g1,A_g2),B_g1),(C_g1,D_g1));")
    species_tree = _mk_tree("((A,B),(C,D));")
    gmap = {
        "A_g1": "A",
        "A_g2": "A",
        "B_g1": "B",
        "C_g1": "C",
        "D_g1": "D",
    }

    res = compute_mulrf_gene_vs_species(gene_tree, species_tree, gene2sp_map=gmap)
    assert res["mulrf_distance"] is not None
    assert res["gene_tree_species_count"] == 4
    assert res["shared_species_count"] == 4


def test_multi_copy_vs_multi_copy_symmetric():
    # Gene-vs-gene MulRF-style conflict should be symmetric.
    tree1 = _mk_tree("(((A_g1,A_g2),B_g1),(C_g1,D_g1));")
    tree2 = _mk_tree("((A_g1,B_g1),(A_g2,(C_g1,D_g1)));")
    gmap = {
        "A_g1": "A",
        "A_g2": "A",
        "B_g1": "B",
        "C_g1": "C",
        "D_g1": "D",
    }

    r12 = compute_mulrf_between_gene_trees(tree1, tree2, gmap, gmap)
    r21 = compute_mulrf_between_gene_trees(tree2, tree1, gmap, gmap)

    assert r12["mulrf_distance"] == r21["mulrf_distance"]
    assert r12["normalized_mulrf_distance"] == r21["normalized_mulrf_distance"]
    assert r12["shared_species_count"] == r21["shared_species_count"]


def test_determinism_same_inputs_repeatable():
    # Re-running should return exactly identical result dictionaries.
    tree1 = _mk_tree("(((A_g1,A_g2),B_g1),(C_g1,D_g1));")
    tree2 = _mk_tree("((A_g1,B_g1),(A_g2,(C_g1,D_g1)));")
    gmap = {
        "A_g1": "A",
        "A_g2": "A",
        "B_g1": "B",
        "C_g1": "C",
        "D_g1": "D",
    }

    first = compute_mulrf_between_gene_trees(tree1, tree2, gmap, gmap)
    for _ in range(5):
        nxt = compute_mulrf_between_gene_trees(tree1, tree2, gmap, gmap)
        assert nxt == first


def test_copy_aware_improves_over_old_species_collapse():
    # Verify that the species-level MulRF distance handles multi-copy genes
    # consistently.  Both tree_a and tree_b have identical species-level
    # bipartitions ({A}, {A,B}, {C,D}), so the copy-aware implementation
    # correctly reports zero distance at the species level.
    tree_a = _mk_tree("((((A_g1,A_g2),A_g3),B_g1),(C_g1,D_g1));")
    tree_b = _mk_tree("((((A_g1,A_g2),B_g1),A_g3),(C_g1,D_g1));")
    gmap = {
        "A_g1": "A",
        "A_g2": "A",
        "A_g3": "A",
        "B_g1": "B",
        "C_g1": "C",
        "D_g1": "D",
    }

    old_d = _old_species_collapse_distance(tree_a, tree_b, gmap, gmap)
    new_res = compute_mulrf_between_gene_trees(tree_a, tree_b, gmap, gmap)

    assert old_d == 0
    assert new_res["mulrf_distance"] is not None
    # Both methods agree: species-level bipartitions are identical.
    assert new_res["mulrf_distance"] == 0


def test_too_few_shared_species_returns_none_distance():
    # Shared species <2 should not crash and should return None distance.
    tree1 = _mk_tree("((A_g1,B_g1));")
    tree2 = _mk_tree("((X_g1,Y_g1));")
    map1 = {"A_g1": "A", "B_g1": "B"}
    map2 = {"X_g1": "X", "Y_g1": "Y"}

    res = compute_mulrf_between_gene_trees(tree1, tree2, map1, map2)
    assert res["shared_species_count"] < 2
    assert res["mulrf_distance"] is None
