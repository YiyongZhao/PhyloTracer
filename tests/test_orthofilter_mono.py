"""
Tests for phylotracer.OrthoFilter_Mono module.
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

from phylotracer.OrthoFilter_Mono import (
    parse_leaf_components,
    get_leaf_clade,
    get_leaf_voucher,
    build_leaf_annotations,
    compute_gene_tree_depths,
    compute_species_tree_depths,
    count_clade_leaves,
    rename_input_single_tre,
)


# =============================================
# parse_leaf_components tests
# =============================================

class TestParseLeafComponents:
    def test_standard_format(self):
        clade, gene_id = parse_leaf_components("Brassica_ABC_1")
        assert clade == "Brassica"
        assert gene_id == "ABC_1"

    def test_two_parts(self):
        clade, gene_id = parse_leaf_components("Clade_Gene1")
        assert clade == "Clade"
        assert gene_id == "Gene1"

    def test_no_underscore(self):
        clade, gene_id = parse_leaf_components("SingleName")
        assert clade == ""
        assert gene_id == "SingleName"

    def test_empty_string(self):
        clade, gene_id = parse_leaf_components("")
        assert clade == ""
        assert gene_id == ""

    def test_multiple_underscores(self):
        clade, gene_id = parse_leaf_components("A_B_C_D")
        assert clade == "A"
        assert gene_id == "B_C_D"


# =============================================
# get_leaf_clade tests
# =============================================

class TestGetLeafClade:
    def test_with_clade_prefix(self):
        result = get_leaf_clade("Brassica_GeneX", {}, {})
        assert result == "Brassica"

    def test_fallback_to_taxa_dic(self):
        result = get_leaf_clade("GeneOnly", {"GeneOnly": "CladeA"}, {})
        assert result == "CladeA"

    def test_fallback_with_renamed_mapping(self):
        result = get_leaf_clade("RenamedGene", {"OrigGene": "CladeB"}, {"RenamedGene": "OrigGene"})
        assert result == "CladeB"

    def test_no_clade_found(self):
        result = get_leaf_clade("Unknown", {}, {})
        assert result == ""

    def test_clade_from_composite_name(self):
        result = get_leaf_clade("MySpecies_Abc_1", {}, {})
        assert result == "MySpecies"


# =============================================
# get_leaf_voucher tests
# =============================================

class TestGetLeafVoucher:
    def test_three_parts(self):
        result = get_leaf_voucher("Clade_VCH_1")
        assert result == "VCH"

    def test_two_parts(self):
        result = get_leaf_voucher("VCH_1")
        assert result == "VCH"

    def test_single_part(self):
        result = get_leaf_voucher("VCH")
        assert result == "VCH"

    def test_many_parts(self):
        result = get_leaf_voucher("A_B_C_D")
        assert result == "C"

    def test_empty_string(self):
        result = get_leaf_voucher("")
        assert result == ""


# =============================================
# build_leaf_annotations tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestBuildLeafAnnotations:
    def test_annotations_returned(self):
        tree = Tree("((X_ABC_1,Y_DEF_2),(Z_GHI_3,W_JKL_4));")
        leaf_to_clade, leaf_to_voucher = build_leaf_annotations(tree, {}, {})
        assert len(leaf_to_clade) == 4
        assert len(leaf_to_voucher) == 4

    def test_clade_values(self):
        tree = Tree("(Sp1_v1_1,Sp2_v2_2);")
        leaf_to_clade, _ = build_leaf_annotations(tree, {}, {})
        assert leaf_to_clade["Sp1_v1_1"] == "Sp1"
        assert leaf_to_clade["Sp2_v2_2"] == "Sp2"


# =============================================
# compute_gene_tree_depths tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestComputeGeneTreeDepths:
    def test_root_depth_zero(self):
        tree = Tree("((A,B),(C,D));")
        depths = compute_gene_tree_depths(tree)
        assert depths[tree] == 0

    def test_leaf_depth(self):
        tree = Tree("((A,B),(C,D));")
        depths = compute_gene_tree_depths(tree)
        leaf = tree & "A"
        assert depths[leaf] == 2

    def test_all_nodes_have_depths(self):
        tree = Tree("((A,B),(C,D));")
        depths = compute_gene_tree_depths(tree)
        for node in tree.traverse():
            assert node in depths


# =============================================
# compute_species_tree_depths tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestComputeSpeciesTreeDepths:
    def test_root_depth_zero(self):
        st = Tree("((A,B),(C,D));")
        depths = compute_species_tree_depths(st)
        assert depths[st] == 0

    def test_consistent_with_gene_tree_depths(self):
        st = Tree("((A,B),(C,D));")
        depths = compute_species_tree_depths(st)
        leaf = st & "A"
        assert depths[leaf] == 2


# =============================================
# count_clade_leaves tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestCountCladeLeaves:
    def test_all_target_clade(self):
        tree = Tree("((A,B),(C,D));")
        leaf_to_clade = {"A": "X", "B": "X", "C": "X", "D": "X"}
        target, total = count_clade_leaves(tree, leaf_to_clade, "X")
        assert target == 4
        assert total == 4

    def test_mixed_clades(self):
        tree = Tree("((A,B),(C,D));")
        leaf_to_clade = {"A": "X", "B": "Y", "C": "X", "D": "Y"}
        target, total = count_clade_leaves(tree, leaf_to_clade, "X")
        assert target == 2
        assert total == 4

    def test_no_target_clade(self):
        tree = Tree("((A,B),(C,D));")
        leaf_to_clade = {"A": "Y", "B": "Y", "C": "Y", "D": "Y"}
        target, total = count_clade_leaves(tree, leaf_to_clade, "X")
        assert target == 0
        assert total == 4

    def test_subtree_count(self):
        tree = Tree("((A,B),(C,D));")
        leaf_to_clade = {"A": "X", "B": "X", "C": "Y", "D": "Y"}
        left_child = tree.get_children()[0]
        target, total = count_clade_leaves(left_child, leaf_to_clade, "X")
        assert target == 2
        assert total == 2


# =============================================
# rename_input_single_tre tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestRenameInputSingleTre:
    def test_renames_leaves(self):
        tree = Tree("(gene1,gene2);")
        taxa_dic = {"gene1": "SpeciesA", "gene2": "SpeciesB"}
        new2gene = {"gene1": "gene1", "gene2": "gene2"}
        result = rename_input_single_tre(tree, taxa_dic, new2gene)
        names = set(result.get_leaf_names())
        assert "SpeciesA_gene1" in names
        assert "SpeciesB_gene2" in names

    def test_missing_taxa_unchanged(self):
        tree = Tree("(gene1,gene2);")
        taxa_dic = {"gene1": "SpeciesA"}
        new2gene = {"gene1": "gene1", "gene2": "gene2"}
        result = rename_input_single_tre(tree, taxa_dic, new2gene)
        names = result.get_leaf_names()
        assert "gene2" in names
