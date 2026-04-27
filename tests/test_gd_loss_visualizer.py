from phylotracer.GD_Loss_Visualizer import (
    get_cumulative_counts_from_paths,
    get_node_event_counts,
    infer_visualizer_mode,
    load_node_summary_counts,
)


def test_get_cumulative_counts_reads_gd_loss_path_column(tmp_path):
    output = tmp_path / "gd_loss.csv"
    output.write_text(
        "\n".join(
            [
                "Tree ID,GD ID,GD burst node,GD loss path",
                'OG_1,7,S20,"S20(2)->S19(2)->S18(1)->S17(0)"',
            ]
        ),
        encoding="utf-8",
    )

    counts = get_cumulative_counts_from_paths(str(output))
    assert counts["S18"]["2_1"] == 1
    assert counts["S17"]["1_0"] == 1


def test_infer_visualizer_mode_prefers_neighbor_node_summary(tmp_path):
    output = tmp_path / "gd_loss.csv"
    output.write_text("Tree ID,GD ID,GD burst node\n", encoding="utf-8")

    node_summary = tmp_path / "gd_loss_node_summary.tsv"
    node_summary.write_text(
        "\n".join(
            [
                "Node\tMode\tGD event count\tNode loss count\tC2_0\tC2_1\tC1_0",
                "S19\taccumulate\t79\t57\t0\t57\t0",
            ]
        ),
        encoding="utf-8",
    )

    assert infer_visualizer_mode(str(output)) == "accumulate"


def test_load_node_summary_counts_reads_accumulate_columns(tmp_path):
    node_summary = tmp_path / "gd_loss_node_summary.tsv"
    node_summary.write_text(
        "\n".join(
            [
                "Node\tMode\tGD event count\tNode loss count\tC2_0\tC2_1\tC1_0",
                "S19\taccumulate\t79\t57\t0\t57\t0",
                "S18\taccumulate\t75\t18\t0\t18\t0",
            ]
        ),
        encoding="utf-8",
    )

    counts, mode, used = load_node_summary_counts(str(node_summary), requested_mode="parsimony")
    assert used is True
    assert mode == "accumulate"
    assert counts["S19"] == {"2_0": 0, "2_1": 57, "1_0": 0, "total": 57}
    assert counts["S18"] == {"2_0": 0, "2_1": 18, "1_0": 0, "total": 18}


def test_get_node_event_counts_uses_gd_loss_pattern_for_parsimony(tmp_path):
    output = tmp_path / "gd_loss.csv"
    output.write_text(
        "\n".join(
            [
                "Tree ID,GD ID,GD burst node,Species,GD loss pattern",
                "OG_1,7,S20,Arabidopsis_thaliana,S18:2_1",
                "OG_1,7,S20,Brassica_rapa,S18:2_1",
                "OG_2,1,S19,Oryza_sativa,S18:2_1;S17:1_0",
            ]
        ),
        encoding="utf-8",
    )

    counts, used = get_node_event_counts(str(output), mode="parsimony")
    assert used is True
    assert counts["S18"]["2_1"] == 2
    assert counts["S17"]["1_0"] == 1


def test_get_node_event_counts_does_not_use_gd_loss_pattern_for_accumulate(tmp_path):
    output = tmp_path / "gd_loss.csv"
    output.write_text(
        "\n".join(
            [
                "Tree ID,GD ID,GD burst node,Species,GD loss pattern",
                "OG_1,7,S20,Arabidopsis_thaliana,S18:2_1",
                "OG_1,7,S20,Brassica_rapa,S18:2_1",
            ]
        ),
        encoding="utf-8",
    )

    counts, used = get_node_event_counts(str(output), mode="accumulate")
    assert used is False
    assert counts == {}
