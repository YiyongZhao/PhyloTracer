"""
Tests for phylotracer.TreeTopology_Summarizer module.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False

from phylotracer.TreeTopology_Summarizer import (
    get_only_sps_tree,
    get_max_tree,
    group_trees_by_topology_with_ids,
    process_tree_with_ids,
)


# =============================================
# get_only_sps_tree tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetOnlySpsTree:
    def test_simplifies_leaf_names(self):
        tree = Tree("((A_1,B_1),(C_1,D_1));")
        result = get_only_sps_tree(tree)
        leaf_names = set(result.get_leaf_names())
        assert leaf_names == {"A", "B", "C", "D"}

    def test_preserves_original_tree(self):
        tree = Tree("((A_1,B_1),(C_1,D_1));")
        get_only_sps_tree(tree)
        assert "A_1" in tree.get_leaf_names()

    def test_handles_multiple_underscores(self):
        tree = Tree("((A_gene_1,B_gene_2),(C_3,D_4));")
        result = get_only_sps_tree(tree)
        leaf_names = set(result.get_leaf_names())
        assert leaf_names == {"A", "B", "C", "D"}

    def test_no_underscore_in_names(self):
        tree = Tree("((Alpha,Beta),(Gamma,Delta));")
        result = get_only_sps_tree(tree)
        leaf_names = set(result.get_leaf_names())
        assert leaf_names == {"Alpha", "Beta", "Gamma", "Delta"}

    def test_single_leaf(self):
        tree = Tree("(X_1);")
        result = get_only_sps_tree(tree)
        assert result.get_leaf_names() == ["X"]

    def test_tree_topology_preserved(self):
        tree = Tree("((A_1,B_1),(C_1,D_1));")
        result = get_only_sps_tree(tree)
        # Should still have two children at root
        assert len(result.get_children()) == 2


# =============================================
# get_max_tree tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetMaxTree:
    def test_returns_largest_tree(self):
        t1 = Tree("(A,B);")
        t2 = Tree("((A,B),(C,D));")
        t3 = Tree("(X,Y,Z);")
        result = get_max_tree([t1, t2, t3])
        assert len(result.get_leaves()) == 4

    def test_single_tree(self):
        t = Tree("(A,B,C);")
        result = get_max_tree([t])
        assert result is t

    def test_empty_list_returns_none(self):
        result = get_max_tree([])
        assert result is None

    def test_equal_size_trees(self):
        t1 = Tree("(A,B);")
        t2 = Tree("(C,D);")
        result = get_max_tree([t1, t2])
        assert len(result.get_leaves()) == 2

    def test_largest_among_many(self):
        trees = [Tree(f"({','.join(chr(65+j) for j in range(i+1))});") for i in range(1, 5)]
        result = get_max_tree(trees)
        assert len(result.get_leaves()) == 5


# =============================================
# group_trees_by_topology_with_ids tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGroupTreesByTopology:
    def test_identical_topologies_grouped(self):
        t1 = Tree("((A,B),(C,D));")
        t2 = Tree("((A,B),(C,D));")
        result_dict = {}
        remaining = group_trees_by_topology_with_ids(
            [("id1", t1), ("id2", t2)], result_dict
        )
        assert len(remaining) == 0
        assert len(result_dict) == 1
        key = list(result_dict.keys())[0]
        assert len(result_dict[key]) == 2

    def test_different_topologies_separated(self):
        t1 = Tree("((A,B),(C,D));")
        t2 = Tree("((A,C),(B,D));")
        result_dict = {}
        remaining = group_trees_by_topology_with_ids(
            [("id1", t1), ("id2", t2)], result_dict
        )
        assert len(remaining) == 1
        assert len(result_dict) == 1

    def test_single_tree_creates_group(self):
        t1 = Tree("((A,B),(C,D));")
        result_dict = {}
        remaining = group_trees_by_topology_with_ids(
            [("id1", t1)], result_dict
        )
        assert len(remaining) == 0
        assert len(result_dict) == 1


# =============================================
# process_tree_with_ids tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestProcessTreeWithIds:
    def test_all_identical_trees(self):
        trees = [("id1", Tree("((A,B),(C,D));")),
                 ("id2", Tree("((A,B),(C,D));")),
                 ("id3", Tree("((A,B),(C,D));"))]
        result_dict = {}
        process_tree_with_ids(trees, result_dict)
        assert len(result_dict) == 1

    def test_all_different_trees(self):
        trees = [("id1", Tree("((A,B),(C,D));")),
                 ("id2", Tree("((A,C),(B,D));")),
                 ("id3", Tree("((A,D),(B,C));"))]
        result_dict = {}
        process_tree_with_ids(trees, result_dict)
        assert len(result_dict) == 3

    def test_mixed_topologies(self):
        trees = [("id1", Tree("((A,B),(C,D));")),
                 ("id2", Tree("((A,B),(C,D));")),
                 ("id3", Tree("((A,C),(B,D));"))]
        result_dict = {}
        process_tree_with_ids(trees, result_dict)
        assert len(result_dict) == 2
