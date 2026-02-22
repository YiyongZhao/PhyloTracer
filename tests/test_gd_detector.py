"""
Tests for phylotracer.GD_Detector module.
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

from phylotracer.GD_Detector import normalize_model, get_model, order_children_by_name

# Import shared utilities needed for annotation
from phylotracer import (
    annotate_gene_tree,
    num_tre_node,
    get_species_set,
    find_dup_node,
)


# =============================================
# normalize_model tests
# =============================================

class TestNormalizeModel:
    def test_aabb(self):
        assert normalize_model("AABB") == "AABB"

    def test_axbb(self):
        assert normalize_model("AXBB") == "AXBB"

    def test_xabb(self):
        assert normalize_model("XABB") == "AXBB"

    def test_aabx(self):
        assert normalize_model("AABX") == "AABX"

    def test_aaxb(self):
        assert normalize_model("AAXB") == "AABX"

    def test_complex_xxxx(self):
        assert normalize_model("XXXX") == "Complex"

    def test_complex_abab(self):
        assert normalize_model("ABAB") == "Complex"

    def test_complex_random(self):
        assert normalize_model("XBXA") == "Complex"


# =============================================
# order_children_by_name tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestOrderChildrenByName:
    def test_orders_by_num_feature(self):
        tree = Tree("((A,B),(C,D));")
        num_tre_node(tree)
        root = tree
        c1, c2 = order_children_by_name(root)
        assert int(c1.num) <= int(c2.num)

    def test_consistent_ordering(self):
        tree = Tree("((A,B),(C,D));")
        num_tre_node(tree)
        c1a, c2a = order_children_by_name(tree)
        c1b, c2b = order_children_by_name(tree)
        assert c1a is c1b
        assert c2a is c2b


# =============================================
# get_model tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestGetModel:
    def test_aabb_pattern(self):
        """Both child clades contain species from both sides of species tree."""
        gene_tree = PhyloTree("((A_1:1,B_1:1)0.9:1,(A_2:1,B_2:1)0.9:1)1.0;", format=1)
        sp_tree = PhyloTree("(A:1,B:1);")
        num_tre_node(gene_tree)
        num_tre_node(sp_tree)
        annotate_gene_tree(gene_tree, sp_tree)
        model = get_model(gene_tree, sp_tree)
        assert model == "AABB"

    def test_axbb_pattern(self):
        """Left clade has only A-species, right has both."""
        gene_tree = PhyloTree("((A_1:1,A_2:1)0.9:1,(A_3:1,B_1:1)0.9:1)1.0;", format=1)
        sp_tree = PhyloTree("(A:1,B:1);")
        num_tre_node(gene_tree)
        num_tre_node(sp_tree)
        annotate_gene_tree(gene_tree, sp_tree)
        model = get_model(gene_tree, sp_tree)
        # One side has A only, other side has A and B
        assert "A" in model
        assert "B" in model

    def test_model_with_four_species(self):
        gene_tree = PhyloTree("((A_1:1,B_1:1)0.9:1,(C_1:1,D_1:1)0.9:1)1.0;", format=1)
        sp_tree = PhyloTree("((A:1,B:1):1,(C:1,D:1):1);")
        num_tre_node(gene_tree)
        num_tre_node(sp_tree)
        annotate_gene_tree(gene_tree, sp_tree)
        model = get_model(gene_tree, sp_tree)
        assert isinstance(model, str)
        assert len(model) == 4


# =============================================
# Integration: annotate + find_dup_node tests
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestFindDupNodeIntegration:
    def test_detects_duplication(self):
        """A tree with clear species overlap should have duplication nodes."""
        gene_tree = PhyloTree("((A_1:1,B_1:1)80:1,(A_2:1,C_1:1)80:1)90;", format=0)
        sp_tree = PhyloTree("((A:1,B:1):1,C:1);")
        num_tre_node(gene_tree)
        num_tre_node(sp_tree)
        annotate_gene_tree(gene_tree, sp_tree)
        dup_nodes = find_dup_node(gene_tree, sp_tree, gd_support=50, clade_support=50,
                                  max_topology_distance=2)
        # Root node should be detected as duplication (A in both children)
        assert len(dup_nodes) >= 1

    def test_no_duplication_in_speciation_tree(self):
        """A tree with no species overlap should have no duplication nodes."""
        gene_tree = PhyloTree("((A_1:1,B_1:1)80:1,(C_1:1,D_1:1)80:1)90;", format=0)
        sp_tree = PhyloTree("((A:1,B:1):1,(C:1,D:1):1);")
        num_tre_node(gene_tree)
        num_tre_node(sp_tree)
        annotate_gene_tree(gene_tree, sp_tree)
        dup_nodes = find_dup_node(gene_tree, sp_tree, gd_support=50, clade_support=50)
        assert len(dup_nodes) == 0

    def test_low_support_not_detected(self):
        """A duplication node with low support should be filtered out."""
        gene_tree = PhyloTree("((A_1:1,B_1:1)10:1,(A_2:1,C_1:1)10:1)20;", format=0)
        sp_tree = PhyloTree("((A:1,B:1):1,C:1);")
        num_tre_node(gene_tree)
        num_tre_node(sp_tree)
        annotate_gene_tree(gene_tree, sp_tree)
        dup_nodes = find_dup_node(gene_tree, sp_tree, gd_support=50, clade_support=50)
        assert len(dup_nodes) == 0

    def test_dup_species_num_filter(self):
        """Requiring more overlap species filters out small duplications."""
        gene_tree = PhyloTree("((A_1:1,B_1:1)80:1,(A_2:1,C_1:1)80:1)90;", format=0)
        sp_tree = PhyloTree("((A:1,B:1):1,C:1);")
        num_tre_node(gene_tree)
        num_tre_node(sp_tree)
        annotate_gene_tree(gene_tree, sp_tree)
        # Require at least 3 overlap species - only A overlaps, so should be filtered
        dup_nodes = find_dup_node(
            gene_tree, sp_tree, gd_support=50, clade_support=50,
            dup_species_num=3
        )
        assert len(dup_nodes) == 0
