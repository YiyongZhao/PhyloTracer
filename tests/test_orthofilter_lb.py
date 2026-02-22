"""
Tests for phylotracer.OrthoFilter_LB module.
"""

import io
import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False

from phylotracer.OrthoFilter_LB import (
    has_multiple_copies,
    get_average_tip_length,
    get_average_node_length,
    remove_long_branches,
)


# =============================================
# has_multiple_copies tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestHasMultipleCopies:
    def test_no_duplicates(self):
        tree = Tree("((A_1,B_1),(C_1,D_1));")
        assert has_multiple_copies(tree) is False

    def test_has_duplicates(self):
        tree = Tree("((A_1,A_2),(B_1,C_1));")
        assert has_multiple_copies(tree) is True

    def test_all_same_species(self):
        tree = Tree("(A_1,A_2,A_3);")
        assert has_multiple_copies(tree) is True

    def test_single_leaf(self):
        tree = Tree("A_1;")
        assert has_multiple_copies(tree) is False


# =============================================
# get_average_tip_length tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetAverageTipLength:
    def test_uniform_branch_lengths(self):
        tree = Tree("((A:2,B:2):1,(C:2,D:2):1);")
        avg = get_average_tip_length(tree)
        assert avg == pytest.approx(2.0)

    def test_varied_branch_lengths(self):
        tree = Tree("((A:1,B:3):1,(C:2,D:4):1);")
        avg = get_average_tip_length(tree)
        assert avg == pytest.approx(2.5)

    def test_zero_branch_lengths(self):
        tree = Tree("((A:0,B:0):0,(C:0,D:0):0);")
        avg = get_average_tip_length(tree)
        assert avg == pytest.approx(0.0)

    def test_single_leaf(self):
        tree = Tree("A:5;")
        avg = get_average_tip_length(tree)
        assert avg == pytest.approx(5.0)


# =============================================
# get_average_node_length tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetAverageNodeLength:
    def test_subtree_average(self):
        tree = Tree("((A:1,B:1):2,(C:3,D:3):2);")
        left_child = tree.get_children()[0]
        avg = get_average_node_length(left_child)
        # distance from subtree root to A = 1, plus subtree.dist = 2
        # distance from subtree root to B = 1, plus subtree.dist = 2
        # total = (1+2 + 1+2) / 2 = 3.0
        assert avg == pytest.approx(3.0)

    def test_leaf_node_average(self):
        tree = Tree("((A:1,B:1):2,(C:3,D:3):2);")
        leaf = tree & "A"
        # For a leaf, subtree has one leaf (itself), dist=1
        # distance from leaf to leaf = 0, plus leaf.dist = 1
        avg = get_average_node_length(leaf)
        assert avg == pytest.approx(1.0)


# =============================================
# remove_long_branches tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestRemoveLongBranches:
    def _make_log(self):
        return io.StringIO()

    def test_no_removal_when_branches_normal(self):
        tree = Tree("((A_1:1,B_1:1):1,(C_1:1,D_1:1):1);")
        log = self._make_log()
        result = remove_long_branches(tree, 5, 2.5, log, "test_tree", {})
        assert set(result.get_leaf_names()) == {"A_1", "B_1", "C_1", "D_1"}

    def test_removes_very_long_branch(self):
        tree = Tree("((A_1:1,B_1:100):1,(C_1:1,D_1:1):1);")
        log = self._make_log()
        result = remove_long_branches(tree, 2, 1, log, "test_tree", {})
        # B_1 should be removed due to extremely long branch
        assert "B_1" not in result.get_leaf_names()

    def test_preserves_original_tree(self):
        tree = Tree("((A_1:1,B_1:100):1,(C_1:1,D_1:1):1);")
        log = self._make_log()
        remove_long_branches(tree, 2, 1, log, "test_tree", {})
        assert "B_1" in tree.get_leaf_names()

    def test_log_records_pruning(self):
        tree = Tree("((A_1:1,B_1:100):1,(C_1:1,D_1:1):1);")
        log = self._make_log()
        remove_long_branches(tree, 2, 1, log, "tree1", {"B_1": "Gene_B"})
        log_content = log.getvalue()
        assert "tree1" in log_content

    def test_zero_dist_leaf_not_removed(self):
        tree = Tree("((A_1:0,B_1:1):1,(C_1:1,D_1:1):1);")
        log = self._make_log()
        result = remove_long_branches(tree, 5, 2.5, log, "test", {})
        assert "A_1" in result.get_leaf_names()

    def test_uses_gene_name_mapping_in_log(self):
        tree = Tree("((A_1:1,B_1:1):1,(C_1:1,D_1:1):1);")
        log = self._make_log()
        mapping = {"A_1": "OrigA", "B_1": "OrigB"}
        remove_long_branches(tree, 5, 2.5, log, "t1", mapping)
        log_content = log.getvalue()
        assert "OrigA" in log_content or "OrigB" in log_content

    def test_all_equal_no_removal(self):
        tree = Tree("((A_1:2,B_1:2):1,(C_1:2,D_1:2):1);")
        log = self._make_log()
        result = remove_long_branches(tree, 5, 2.5, log, "t1", {})
        assert len(result.get_leaves()) == 4
