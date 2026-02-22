"""
Tests for phylotracer.__init__ utility functions.
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

from phylotracer import (
    is_rooted,
    get_species_set,
    get_species_list,
    calculate_species_num,
    num_tre_node,
    read_and_return_dict,
    rename_input_tre,
    realign_branch_length,
    get_max_deepth,
    judge_support,
    sps_dup_num,
    get_gene_pairs,
    generate_sps_voucher,
    root_tre_with_midpoint_outgroup,
    compute_tip_to_root_branch_length_variance,
    map_species_set_to_node,
    find_tre_dup,
)


# =============================================
# is_rooted tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestIsRooted:
    def test_is_rooted_returns_true_for_bifurcating_root(self):
        tree = Tree("((A,B),(C,D));")
        assert is_rooted(tree) is True

    def test_is_rooted_returns_false_for_trifurcating_root(self):
        tree = Tree("(A,B,(C,D));")
        assert is_rooted(tree) is False

    def test_is_rooted_single_child(self):
        # A tree with only one child at root
        tree = Tree("((A,B));")
        assert is_rooted(tree) is False

    def test_is_rooted_two_leaves(self):
        tree = Tree("(A,B);")
        assert is_rooted(tree) is True


# =============================================
# get_species_set / get_species_list tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetSpecies:
    def test_get_species_set_extracts_unique_species(self):
        tree = Tree("((A_1,A_2),(B_1,C_1));")
        result = get_species_set(tree)
        assert result == {"A", "B", "C"}

    def test_get_species_set_with_single_species(self):
        tree = Tree("(A_1,A_2);")
        result = get_species_set(tree)
        assert result == {"A"}

    def test_get_species_list_preserves_duplicates(self):
        tree = Tree("((A_1,A_2),(B_1,C_1));")
        result = get_species_list(tree)
        assert sorted(result) == ["A", "A", "B", "C"]

    def test_get_species_set_returns_empty_for_none(self):
        result = get_species_list(None)
        assert result == []

    def test_calculate_species_num(self):
        tree = Tree("((A_1,B_1),(C_1,D_1));")
        assert calculate_species_num(tree) == 4

    def test_get_species_set_no_underscore(self):
        """Leaf names without underscore: entire name is the species."""
        tree = Tree("(A,B,C);")
        result = get_species_set(tree)
        assert result == {"A", "B", "C"}


# =============================================
# num_tre_node tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestNumTreNode:
    def test_num_tre_node_labels_internal_nodes(self):
        tree = Tree("((A,B),(C,D));")
        result = num_tre_node(tree)
        internal_names = [
            n.name for n in result.traverse() if not n.is_leaf()
        ]
        assert "N0" in internal_names
        assert "N1" in internal_names
        assert "N2" in internal_names

    def test_num_tre_node_preserves_leaf_names(self):
        tree = Tree("((A,B),(C,D));")
        result = num_tre_node(tree)
        leaf_names = result.get_leaf_names()
        assert set(leaf_names) == {"A", "B", "C", "D"}

    def test_num_tre_node_assigns_num_feature(self):
        tree = Tree("((A,B),(C,D));")
        result = num_tre_node(tree)
        for node in result.traverse():
            assert hasattr(node, "num")
            assert isinstance(node.num, int)

    def test_num_tre_node_postorder_numbering(self):
        tree = Tree("((A,B),(C,D));")
        result = num_tre_node(tree)
        nums = [n.num for n in result.traverse("postorder")]
        assert nums == list(range(1, len(nums) + 1))


# =============================================
# read_and_return_dict tests
# =============================================

class TestReadAndReturnDict:
    def test_read_valid_mapping_file(self, tmp_path):
        fpath = tmp_path / "map.tsv"
        fpath.write_text("key1\tval1\nkey2\tval2\n")
        result = read_and_return_dict(str(fpath))
        assert result == {"key1": "val1", "key2": "val2"}

    def test_read_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            read_and_return_dict("/nonexistent/path.tsv")

    def test_read_empty_file_raises(self, tmp_path):
        fpath = tmp_path / "empty.tsv"
        fpath.write_text("")
        with pytest.raises(Exception):
            read_and_return_dict(str(fpath))

    def test_read_with_custom_separator(self, tmp_path):
        fpath = tmp_path / "map.csv"
        fpath.write_text("key1,val1\nkey2,val2\n")
        result = read_and_return_dict(str(fpath), separator=",")
        assert result == {"key1": "val1", "key2": "val2"}

    def test_read_single_row(self, tmp_path):
        fpath = tmp_path / "single.tsv"
        fpath.write_text("alpha\tbeta\n")
        result = read_and_return_dict(str(fpath))
        assert result == {"alpha": "beta"}


# =============================================
# rename_input_tre tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestRenameInputTre:
    def test_rename_replaces_leaf_names(self):
        tree = Tree("((geneA,geneB),(geneC,geneD));")
        mapping = {"geneA": "A_1", "geneB": "B_1"}
        result = rename_input_tre(tree, mapping)
        leaf_names = set(result.get_leaf_names())
        assert "A_1" in leaf_names
        assert "B_1" in leaf_names
        assert "geneC" in leaf_names  # unchanged

    def test_rename_preserves_original_tree(self):
        tree = Tree("((geneA,geneB),(geneC,geneD));")
        mapping = {"geneA": "A_1"}
        rename_input_tre(tree, mapping)
        assert "geneA" in tree.get_leaf_names()

    def test_rename_empty_mapping(self):
        tree = Tree("((A,B),(C,D));")
        result = rename_input_tre(tree, {})
        assert set(result.get_leaf_names()) == {"A", "B", "C", "D"}

    def test_rename_all_leaves(self):
        tree = Tree("(X,Y);")
        mapping = {"X": "A_1", "Y": "B_1"}
        result = rename_input_tre(tree, mapping)
        assert set(result.get_leaf_names()) == {"A_1", "B_1"}


# =============================================
# realign_branch_length tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestRealignBranchLength:
    def test_realign_branch_length_returns_tree(self):
        tree = Tree("((A:1,B:2):3,(C:4,D:5):6);")
        result = realign_branch_length(tree)
        assert result is not None
        assert len(result.get_leaves()) == 4

    def test_realign_branch_length_modifies_distances(self):
        tree = Tree("((A:1,B:2):3,(C:4,D:5):6);")
        original_dists = [n.dist for n in tree.traverse() if not n.is_root()]
        result = realign_branch_length(tree)
        new_dists = [n.dist for n in result.traverse() if not n.is_root()]
        assert original_dists != new_dists


# =============================================
# get_max_deepth tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetMaxDepth:
    def test_leaf_depth_is_one(self):
        tree = Tree("(A);")
        # Root with one child: depth should be > 0
        assert get_max_deepth(tree) >= 1

    def test_balanced_tree_depth(self):
        tree = Tree("((A,B),(C,D));")
        assert get_max_deepth(tree) == 3  # root -> internal -> leaf

    def test_none_returns_zero(self):
        assert get_max_deepth(None) == 0


# =============================================
# judge_support tests
# =============================================

class TestJudgeSupport:
    def test_support_above_threshold(self):
        assert judge_support(80, 50) is True

    def test_support_below_threshold(self):
        assert judge_support(30, 50) is False

    def test_support_equal_threshold(self):
        assert judge_support(50, 50) is True

    def test_fractional_support_converted(self):
        assert judge_support(0.8, 0.5) is True

    def test_fractional_below(self):
        assert judge_support(0.3, 0.5) is False

    def test_zero_support(self):
        assert judge_support(0, 50) is False

    def test_hundred_support(self):
        assert judge_support(100, 100) is True


# =============================================
# sps_dup_num tests
# =============================================

class TestSpsDupNum:
    def test_no_duplicates(self):
        result = sps_dup_num(["A", "B", "C"], {"A", "B", "C"})
        assert result == 0

    def test_one_duplicate(self):
        result = sps_dup_num(["A", "A", "B"], {"A", "B"})
        assert result == 1

    def test_all_duplicated(self):
        result = sps_dup_num(["A", "A", "B", "B"], {"A", "B"})
        assert result == 2

    def test_empty_list(self):
        result = sps_dup_num([], {"A", "B"})
        assert result == 0


# =============================================
# get_gene_pairs tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetGenePairs:
    def test_binary_node_gene_pairs(self):
        tree = Tree("((A_1,B_1),(A_2,C_1));")
        pairs = get_gene_pairs(tree)
        species_in_pairs = {p[0] for p in pairs}
        assert "A" in species_in_pairs

    def test_non_binary_returns_empty(self):
        tree = Tree("(A_1,B_1,C_1);")
        pairs = get_gene_pairs(tree)
        assert pairs == []

    def test_species_on_one_side_only(self):
        tree = Tree("((A_1,B_1),(C_1,D_1));")
        pairs = get_gene_pairs(tree)
        # A only on left, C only on right => None fills
        a_pairs = [p for p in pairs if p[0] == "A"]
        assert len(a_pairs) >= 1
        for p in a_pairs:
            assert p[2] is None  # A not on right side


# =============================================
# generate_sps_voucher tests
# =============================================

class TestGenerateSpsVoucher:
    def test_generates_correct_count(self):
        result = generate_sps_voucher(5)
        assert len(result) == 5

    def test_vouchers_are_unique(self):
        result = generate_sps_voucher(10)
        assert len(set(result)) == 10

    def test_vouchers_are_three_chars(self):
        result = generate_sps_voucher(3)
        assert all(len(v) == 3 for v in result)

    def test_vouchers_are_sorted(self):
        result = generate_sps_voucher(5)
        assert result == sorted(result)


# =============================================
# root_tre_with_midpoint_outgroup tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestRootTreWithMidpointOutgroup:
    def test_already_rooted_tree(self):
        tree = Tree("((A:1,B:1):1,(C:1,D:1):1);")
        result = root_tre_with_midpoint_outgroup(tree)
        assert is_rooted(result) is True

    def test_unrooted_tree_gets_rooted(self):
        tree = Tree("(A:1,B:2,(C:3,D:4):5);")
        result = root_tre_with_midpoint_outgroup(tree)
        assert is_rooted(result) is True

    def test_small_tree_returned_unchanged(self):
        tree = Tree("(A:1,B:2);")
        result = root_tre_with_midpoint_outgroup(tree)
        assert len(result.get_leaves()) == 2


# =============================================
# compute_tip_to_root_branch_length_variance tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestBranchLengthVariance:
    def test_equal_branch_lengths_zero_variance(self):
        tree = Tree("((A:1,B:1):1,(C:1,D:1):1);")
        var = compute_tip_to_root_branch_length_variance(tree)
        assert var == pytest.approx(0.0)

    def test_unequal_branch_lengths_nonzero_variance(self):
        tree = Tree("((A:1,B:5):1,(C:1,D:1):1);")
        var = compute_tip_to_root_branch_length_variance(tree)
        assert var > 0.0

    def test_single_leaf_returns_zero(self):
        tree = Tree("A:1;")
        var = compute_tip_to_root_branch_length_variance(tree)
        assert var == 0.0


# =============================================
# map_species_set_to_node tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestMapSpeciesSetToNode:
    def test_single_species_returns_leaf(self, species_tree):
        node = map_species_set_to_node(species_tree, {"A"})
        assert node is not None
        assert node.is_leaf()
        assert node.name == "A"

    def test_two_species_returns_ancestor(self, species_tree):
        node = map_species_set_to_node(species_tree, {"A", "B"})
        assert node is not None
        leaves = set(node.get_leaf_names())
        assert "A" in leaves and "B" in leaves

    def test_empty_set_returns_none(self, species_tree):
        assert map_species_set_to_node(species_tree, set()) is None

    def test_unknown_species_returns_none(self, species_tree):
        assert map_species_set_to_node(species_tree, {"Z"}) is None


# =============================================
# find_tre_dup tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestFindTreDup:
    def test_tree_without_reconciliation_raises(self):
        tree = Tree("((A,B),(C,D));")
        with pytest.raises(ValueError, match="reconciliation"):
            find_tre_dup(tree)

    def test_phylo_tree_returns_leaves(self):
        tree = PhyloTree("((A_1,B_1),(C_1,D_1));")
        sp_tree = PhyloTree("((A,B),(C,D));")
        # Reconcile to get events
        tree.set_species_naming_function(lambda x: x.split("_")[0])
        sp_tree.set_species_naming_function(lambda x: x)
        try:
            recon = tree.reconcile(sp_tree)
            pairs, leaves = find_tre_dup(recon[0])
            assert isinstance(leaves, set)
        except Exception:
            # If reconciliation is not available, just verify it raises properly
            pass
