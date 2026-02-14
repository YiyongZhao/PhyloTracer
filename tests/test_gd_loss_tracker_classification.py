import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))

from GD_Loss_Tracker import classify_species_copy_state


class FakeLeaf:
    def __init__(self, name):
        self.name = name
        self.children = []

    def iter_leaves(self):
        yield self

    def get_leaf_names(self):
        return [self.name]


class FakeNode:
    def __init__(self, name, children=None):
        self.name = name
        self.children = children or []

    def get_children(self):
        return self.children

    def iter_leaves(self):
        if not self.children:
            yield FakeLeaf(self.name)
            return
        for child in self.children:
            yield from child.iter_leaves()

    def get_leaf_names(self):
        return [leaf.name for leaf in self.iter_leaves()]


def make_binary_dup(left_leaf_names, right_leaf_names):
    left = FakeNode("L", [FakeNode(n) for n in left_leaf_names])
    right = FakeNode("R", [FakeNode(n) for n in right_leaf_names])
    return FakeNode("D", [left, right])


def make_mapped_node(clade_species):
    return FakeNode("Sx", [FakeNode(sp) for sp in clade_species])


def test_case_a_2_2_both_sides_present():
    dup = make_binary_dup(["A_g1", "B_g1"], ["A_g2", "C_g1"])
    mapped = make_mapped_node(["A", "B", "C"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"A", "B", "C"}
    )
    assert (loss_type, left_has, right_has, confidence) == ("2-2", True, True, "observable")


def test_case_b_2_1_left_only():
    dup = make_binary_dup(["A_g1", "B_g1"], ["C_g1", "D_g1"])
    mapped = make_mapped_node(["A", "B", "C", "D"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"A", "B", "C", "D"}
    )
    assert (loss_type, left_has, right_has, confidence) == ("2-1", True, False, "observable")


def test_case_c_2_0_when_species_in_family_but_absent_in_dup_children():
    dup = make_binary_dup(["B_g1"], ["C_g1"])
    mapped = make_mapped_node(["A", "B", "C"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"A", "B", "C"}
    )
    assert (loss_type, left_has, right_has, confidence) == ("2-0", False, False, "observable")


def test_case_d_missing_data_default_excluded():
    dup = make_binary_dup(["B_g1"], ["C_g1"])
    mapped = make_mapped_node(["A", "B", "C"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"B", "C"}, include_unobserved_species=False
    )
    assert (loss_type, left_has, right_has, confidence) == ("NA", False, False, "unobserved_in_family")


def test_case_e_leaf_label_without_underscore_is_parsed():
    dup = make_binary_dup(["A"], ["A_g2", "B_g1"])
    mapped = make_mapped_node(["A", "B"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"A", "B"}
    )
    assert (loss_type, left_has, right_has, confidence) == ("2-2", True, True, "observable")


def test_case_f_invalid_node_or_mapping_returns_na():
    non_binary_dup = FakeNode("D", [FakeNode("A_g1"), FakeNode("A_g2"), FakeNode("B_g1")])
    mapped = make_mapped_node(["A", "B"])
    res1 = classify_species_copy_state("A", non_binary_dup, mapped, {"A", "B"})
    res2 = classify_species_copy_state("A", make_binary_dup(["A_g1"], ["B_g1"]), None, {"A", "B"})
    assert res1 == ("NA", False, False, "invalid_node")
    assert res2 == ("NA", False, False, "invalid_node")
