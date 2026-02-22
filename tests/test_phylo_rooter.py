"""
Tests for phylotracer.Phylo_Rooter module.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    from ete3 import Tree, PhyloTree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False

from phylotracer.Phylo_Rooter import (
    get_species_map_and_depth,
    get_species_tree_basal_set,
    get_dynamic_basal_set,
    annotate_mapped_depths,
)


# =============================================
# get_species_map_and_depth tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetSpeciesMapAndDepth:
    def test_root_has_depth_zero(self):
        st = Tree("((A,B),(C,D));")
        depths = get_species_map_and_depth(st)
        assert depths[st] == 0

    def test_leaves_have_depth_two(self):
        st = Tree("((A,B),(C,D));")
        depths = get_species_map_and_depth(st)
        leaf_a = st & "A"
        assert depths[leaf_a] == 2

    def test_internal_node_depth(self):
        st = Tree("((A,B),(C,D));")
        depths = get_species_map_and_depth(st)
        # Internal nodes at depth 1
        for child in st.get_children():
            assert depths[child] == 1

    def test_all_nodes_have_depth(self):
        st = Tree("((A,B),(C,D));")
        depths = get_species_map_and_depth(st)
        for node in st.traverse():
            assert node in depths

    def test_asymmetric_tree(self):
        st = Tree("(A,(B,(C,D)));")
        depths = get_species_map_and_depth(st)
        leaf_d = st & "D"
        assert depths[leaf_d] == 3


# =============================================
# get_species_tree_basal_set tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetSpeciesTreeBasalSet:
    def test_returns_smaller_clade(self):
        st = Tree("((A,B,C),(D,E));")
        basal = get_species_tree_basal_set(st)
        assert basal == {"D", "E"}

    def test_balanced_tree(self):
        st = Tree("((A,B),(C,D));")
        basal = get_species_tree_basal_set(st)
        # Both clades have 2 leaves; either is valid
        assert len(basal) == 2

    def test_single_child_returns_all(self):
        st = Tree("(A);")
        basal = get_species_tree_basal_set(st)
        assert "A" in basal

    def test_asymmetric_returns_singleton(self):
        st = Tree("(A,(B,(C,D)));")
        basal = get_species_tree_basal_set(st)
        assert basal == {"A"}


# =============================================
# get_dynamic_basal_set tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetDynamicBasalSet:
    def test_all_species_present(self):
        st = Tree("((A,B),(C,D));")
        gene_sps = {"A", "B", "C", "D"}
        basal = get_dynamic_basal_set(gene_sps, st)
        # Should return one of the root clades
        assert len(basal) > 0
        assert basal.issubset(gene_sps)

    def test_species_only_in_one_clade(self):
        st = Tree("((A,B),(C,D));")
        gene_sps = {"A", "B"}  # only in left clade
        basal = get_dynamic_basal_set(gene_sps, st)
        # Should recurse into the left clade
        assert basal.issubset({"A", "B"})

    def test_empty_species_set(self):
        st = Tree("((A,B),(C,D));")
        basal = get_dynamic_basal_set(set(), st)
        assert basal == set()

    def test_leaf_species_tree(self):
        st = Tree("A;")
        basal = get_dynamic_basal_set({"A"}, st)
        assert basal == {"A"}


# =============================================
# annotate_mapped_depths tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestAnnotateMappedDepths:
    def test_leaves_get_mapped_depth(self):
        gt = Tree("((A_1,B_1),(C_1,D_1));")
        st = Tree("((A,B),(C,D));")
        result = annotate_mapped_depths(gt, st)
        for leaf in result.iter_leaves():
            assert hasattr(leaf, "mapped_depth")

    def test_none_gene_tree_returned(self):
        st = Tree("((A,B),(C,D));")
        assert annotate_mapped_depths(None, st) is None

    def test_none_species_tree_returned(self):
        gt = Tree("((A_1,B_1),(C_1,D_1));")
        result = annotate_mapped_depths(gt, None)
        assert result is gt

    def test_internal_nodes_get_mapped_depth(self):
        gt = Tree("((A_1,B_1),(C_1,D_1));")
        st = Tree("((A,B),(C,D));")
        result = annotate_mapped_depths(gt, st)
        for node in result.traverse():
            assert hasattr(node, "mapped_depth")

    def test_mapped_depth_root_is_zero(self):
        gt = Tree("((A_1,B_1),(C_1,D_1));")
        st = Tree("((A,B),(C,D));")
        result = annotate_mapped_depths(gt, st)
        # Root maps to root of species tree => depth 0
        assert result.mapped_depth == 0
