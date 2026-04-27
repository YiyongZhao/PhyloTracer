import pandas as pd
import pytest

from phylotracer.GD_Loss_Tracker import (
    LOSS_TYPE_2_0,
    LOSS_TYPE_2_1,
    LOSS_TYPE_2_2,
    LOSS_TYPE_NA,
    classify_species_copy_state,
    infer_loss_node_from_path,
    write_gd_loss_csv,
    write_gd_loss_node_summary_tsv,
)


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
    assert (loss_type, left_has, right_has, confidence) == (LOSS_TYPE_2_2, True, True, "observable")


def test_case_b_2_1_left_only():
    dup = make_binary_dup(["A_g1", "B_g1"], ["C_g1", "D_g1"])
    mapped = make_mapped_node(["A", "B", "C", "D"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"A", "B", "C", "D"}
    )
    assert (loss_type, left_has, right_has, confidence) == (LOSS_TYPE_2_1, True, False, "observable")


def test_case_c_2_0_when_species_in_family_but_absent_in_dup_children():
    dup = make_binary_dup(["B_g1"], ["C_g1"])
    mapped = make_mapped_node(["A", "B", "C"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"A", "B", "C"}
    )
    assert (loss_type, left_has, right_has, confidence) == (LOSS_TYPE_2_0, False, False, "observable")


def test_case_d_missing_data_default_excluded():
    dup = make_binary_dup(["B_g1"], ["C_g1"])
    mapped = make_mapped_node(["A", "B", "C"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"B", "C"}, include_unobserved_species=False
    )
    assert (loss_type, left_has, right_has, confidence) == (LOSS_TYPE_NA, False, False, "unobserved_in_family")


def test_case_e_leaf_label_without_underscore_is_parsed():
    dup = make_binary_dup(["A"], ["A_g2", "B_g1"])
    mapped = make_mapped_node(["A", "B"])
    loss_type, left_has, right_has, confidence = classify_species_copy_state(
        "A", dup, mapped, {"A", "B"}
    )
    assert (loss_type, left_has, right_has, confidence) == (LOSS_TYPE_2_2, True, True, "observable")


def test_case_f_invalid_node_or_mapping_returns_na():
    non_binary_dup = FakeNode("D", [FakeNode("A_g1"), FakeNode("A_g2"), FakeNode("B_g1")])
    mapped = make_mapped_node(["A", "B"])
    res1 = classify_species_copy_state("A", non_binary_dup, mapped, {"A", "B"})
    res2 = classify_species_copy_state("A", make_binary_dup(["A_g1"], ["B_g1"]), None, {"A", "B"})
    assert res1 == ("NA", False, False, "invalid_node")
    assert res2 == ("NA", False, False, "invalid_node")


def test_infer_loss_node_from_path_uses_first_copy_drop():
    assert infer_loss_node_from_path("S20(2)->S19(2)->S18(1)->S17(1)") == "S18"
    assert infer_loss_node_from_path("S20(2)->S19(1)->S18(1)") == "S19"
    assert infer_loss_node_from_path("S20(2)->S19(2)->S18(0)") == "S18"
    assert infer_loss_node_from_path("S20(2)->S19(2)") == "NA"
    assert infer_loss_node_from_path("NA") == "NA"


def test_write_gd_loss_csv_uses_clear_node_event_columns_without_loss_node(tmp_path):
    output = tmp_path / "gd_loss.csv"
    write_gd_loss_csv(
        [
            {
                "tre_id": "OG_1",
                "gd_id": "7",
                "gd_level_name": "S20",
                "support": "95.0",
                "species_display": "Arabidopsis_thaliana",
                "gene1": "ARA_A1",
                "gene2": "ARA_B1",
                "left_gene_n": 1,
                "right_gene_n": 0,
                "left_has": 1,
                "right_has": 0,
                "loss_type": "2_1",
                "loss_confidence": "observable",
                "major_loss_class": "loss_one_copy",
                "loss_path": "S20(2)->S19(2)->S18(1)->Arabidopsis_thaliana(1)",
                "effective_transition_keys_for_export": [("S18", "2_1")],
                "path_count_node_events": "S18|2_1",
                "path_count_types": "2_1",
                "path_count_2_0": 0,
                "path_count_2_1": 1,
                "path_count_1_0": 0,
            }
        ],
        str(output),
    )
    df = pd.read_csv(output)
    assert list(df.columns[:5]) == ["Tree ID", "GD ID", "GD burst node", "GD loss node", "GD loss pattern"]
    assert "Loss node" not in df.columns
    assert df.loc[0, "GD subclade A genes"] == "ARA_A1"
    assert df.loc[0, "GD subclade B genes"] == "ARA_B1"
    assert df.loc[0, "GD loss node"] == "S18"
    assert df.loc[0, "GD loss pattern"] == "S18:2_1"
    assert df.columns[-1] == "GD loss path"
    assert df.loc[0, "GD loss path"] == "S20(2)->S19(2)->S18(1)->Arabidopsis_thaliana(1)"


def test_write_gd_loss_csv_keeps_species_specific_subclade_genes(tmp_path):
    output = tmp_path / "gd_loss.csv"
    write_gd_loss_csv(
        [
            {
                "tre_id": "OG_1",
                "gd_id": "7",
                "gd_level_name": "S20",
                "support": "95.0",
                "species_display": "Chlamydomonas_reinhardtii",
                "gene1": "CRE_A1",
                "gene2": "CRE_B1",
                "left_gene_n": 1,
                "right_gene_n": 1,
                "loss_type": "2_2",
                "loss_path": "S20(2)->S1(2)->S0(2)->Chlamydomonas_reinhardtii(2)",
                "effective_transition_keys_for_export": [],
            },
            {
                "tre_id": "OG_1",
                "gd_id": "7",
                "gd_level_name": "S20",
                "support": "95.0",
                "species_display": "Actinidia_chinensis",
                "gene1": "ACT_A1",
                "gene2": "ACT_B1,ACT_B2",
                "left_gene_n": 1,
                "right_gene_n": 2,
                "loss_type": "2_2",
                "loss_path": "S20(2)->S19(2)->S18(2)->Actinidia_chinensis(2)",
                "effective_transition_keys_for_export": [],
            },
        ],
        str(output),
    )
    df = pd.read_csv(output)
    cre_row = df[df["Species"] == "Chlamydomonas_reinhardtii"].iloc[0]
    act_row = df[df["Species"] == "Actinidia_chinensis"].iloc[0]
    assert cre_row["GD subclade A genes"] == "CRE_A1"
    assert cre_row["GD subclade B genes"] == "CRE_B1"
    assert cre_row["GD subclade A gene count"] == 1
    assert cre_row["GD subclade B gene count"] == 1
    assert act_row["GD subclade A genes"] == "ACT_A1"
    assert act_row["GD subclade B genes"] == "ACT_B1,ACT_B2"
    assert act_row["GD subclade A gene count"] == 1
    assert act_row["GD subclade B gene count"] == 2


def test_write_gd_loss_csv_uses_event_level_loss_node_and_pattern(tmp_path):
    output = tmp_path / "gd_loss.csv"
    write_gd_loss_csv(
        [
            {
                "tre_id": "OG_1",
                "gd_id": "7",
                "gd_level_name": "S20",
                "support": "95.0",
                "species_display": "Chlamydomonas_reinhardtii",
                "gene1": "CRE_A1",
                "gene2": "CRE_B1",
                "left_gene_n": 1,
                "right_gene_n": 1,
                "loss_type": "2_2",
                "loss_path": "S20(2)->S1(2)->S0(2)->Chlamydomonas_reinhardtii(2)",
                "effective_transition_keys_for_export": [],
            },
            {
                "tre_id": "OG_1",
                "gd_id": "7",
                "gd_level_name": "S20",
                "support": "95.0",
                "species_display": "Actinidia_chinensis",
                "gene1": "NA",
                "gene2": "ACT_B1,ACT_B2",
                "left_gene_n": 0,
                "right_gene_n": 2,
                "loss_type": "2_1",
                "loss_path": "S20(2)->S19(1)->S18(1)->Actinidia_chinensis(1)",
                "effective_transition_keys_for_export": [("S19", "2_1"), ("Actinidia_chinensis", "1_0")],
            },
        ],
        str(output),
    )
    df = pd.read_csv(output)
    assert df["GD loss node"].tolist() == ["S19", "S19"]
    assert df["GD loss pattern"].tolist() == ["S19:2_1", "S19:2_1"]


def test_node_summary_uses_p_columns_for_parsimony(tmp_path):
    output = tmp_path / "gd_loss_node_summary.tsv"
    records = [
        {
            "tre_id": "OG_1",
            "gd_id": "7",
            "gd_level_name": "S20",
            "effective_transition_keys_for_export": [("S18", "2_1")],
        },
        {
            "tre_id": "OG_1",
            "gd_id": "7",
            "gd_level_name": "S20",
            "effective_transition_keys_for_export": [("S18", "2_1")],
        },
    ]
    write_gd_loss_node_summary_tsv(records, str(output), node_count_mode="parsimony")
    df = pd.read_csv(output, sep="\t")
    assert {"P2_0", "P2_1", "P1_0"} <= set(df.columns)
    assert not {"C2_0", "C2_1", "C1_0"} & set(df.columns)
    s18 = df[df["Node"] == "S18"].iloc[0]
    assert s18["P2_1"] == 1
    assert s18["Node loss count"] == 1


def test_node_summary_uses_c_columns_for_accumulate(tmp_path):
    output = tmp_path / "gd_loss_node_summary.tsv"
    records = [
        {
            "tre_id": "OG_1",
            "gd_id": "7",
            "gd_level_name": "S20",
            "effective_transition_keys_for_export": [("S18", "2_1")],
        },
        {
            "tre_id": "OG_1",
            "gd_id": "7",
            "gd_level_name": "S20",
            "effective_transition_keys_for_export": [("S18", "2_1")],
        },
    ]
    write_gd_loss_node_summary_tsv(records, str(output), node_count_mode="accumulate")
    df = pd.read_csv(output, sep="\t")
    assert {"C2_0", "C2_1", "C1_0"} <= set(df.columns)
    assert not {"P2_0", "P2_1", "P1_0"} & set(df.columns)
    s18 = df[df["Node"] == "S18"].iloc[0]
    assert s18["C2_1"] == 2
    assert s18["Node loss count"] == 2
