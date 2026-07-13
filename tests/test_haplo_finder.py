"""
Tests for phylotracer.HaploFinder module.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from phylotracer.HaploFinder import (
    confident_projection_pairs,
    find_conversion_zones_with_ids_to_file,
    infer_projection_color_map,
    read_gff,
    read_lens,
    split_sequences,
)


def _projection_fixture(donor_chrs, target_chrs, pair_counts, background_genes=100):
    """Build compact chromosome projection data with fixed chromosome sizes."""
    dict_gff1 = {}
    dict_gff2 = {}
    for chromosome in donor_chrs:
        for idx in range(background_genes):
            dict_gff1[f"a_{chromosome}_{idx}"] = [
                chromosome, str(idx + 1), str(idx + 1), "+", ".",
            ]
    for chromosome in target_chrs:
        for idx in range(background_genes):
            dict_gff2[f"b_{chromosome}_{idx}"] = [
                chromosome, str(idx + 1), str(idx + 1), "+", ".",
            ]

    raw_colors = {}
    for (chr_a, chr_b), count in pair_counts.items():
        for idx in range(count):
            raw_colors[f"a_{chr_a}_{idx}:b_{chr_b}_{idx}"] = (
                "red" if idx % 2 == 0 else "blue"
            )
    lens_a = [[chromosome, str(background_genes)] for chromosome in donor_chrs]
    lens_b = [[chromosome, str(background_genes)] for chromosome in target_chrs]
    return raw_colors, dict_gff1, dict_gff2, lens_a, lens_b


# =============================================
# read_gff tests
# =============================================

class TestReadGff:
    def test_reads_valid_gff(self, tmp_path):
        content = (
            "chr1\tgene1\t100\t200\t+\tattr1\n"
            "chr2\tgene2\t300\t400\t-\tattr2\n"
        )
        fpath = tmp_path / "test.gff"
        fpath.write_text(content)
        data, data_dict = read_gff(str(fpath))
        assert len(data) == 2
        assert "gene1" in data_dict
        assert "gene2" in data_dict

    def test_gff_data_list_structure(self, tmp_path):
        content = "chr1\tgene1\t100\t200\t+\tattr1\n"
        fpath = tmp_path / "test.gff"
        fpath.write_text(content)
        data, _ = read_gff(str(fpath))
        assert data[0][0] == "chr1"
        assert data[0][1] == "gene1"

    def test_gff_dict_values(self, tmp_path):
        content = "chr1\tgene1\t100\t200\t+\tattr1\n"
        fpath = tmp_path / "test.gff"
        fpath.write_text(content)
        _, data_dict = read_gff(str(fpath))
        assert data_dict["gene1"] == ["chr1", "100", "200", "+", "attr1"]

    def test_empty_gff(self, tmp_path):
        fpath = tmp_path / "empty.gff"
        fpath.write_text("")
        data, data_dict = read_gff(str(fpath))
        assert data == []
        assert data_dict == {}

    def test_multiple_genes_same_chr(self, tmp_path):
        content = (
            "chr1\tgeneA\t100\t200\t+\ta1\n"
            "chr1\tgeneB\t300\t400\t+\ta2\n"
            "chr1\tgeneC\t500\t600\t-\ta3\n"
        )
        fpath = tmp_path / "test.gff"
        fpath.write_text(content)
        data, data_dict = read_gff(str(fpath))
        assert len(data) == 3
        assert len(data_dict) == 3

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            read_gff("/nonexistent/file.gff")


# =============================================
# read_lens tests
# =============================================

class TestReadLens:
    def test_reads_valid_lens(self, tmp_path):
        content = "chr1\t1000\nchr2\t2000\nchr3\t3000\n"
        fpath = tmp_path / "test.lens"
        fpath.write_text(content)
        data = read_lens(str(fpath))
        assert len(data) == 3
        assert data[0] == ["chr1", "1000"]
        assert data[1] == ["chr2", "2000"]

    def test_lens_with_chr_filter(self, tmp_path):
        lens_content = "chr1\t1000\nchr2\t2000\nchr3\t3000\n"
        fpath = tmp_path / "test.lens"
        fpath.write_text(lens_content)

        chrs_content = "chr1\nchr3\n"
        chrs_path = tmp_path / "chrs.txt"
        chrs_path.write_text(chrs_content)

        data = read_lens(str(fpath), chrs=str(chrs_path))
        assert len(data) == 2
        chr_names = [d[0] for d in data]
        assert "chr1" in chr_names
        assert "chr3" in chr_names
        assert "chr2" not in chr_names

    def test_lens_no_filter(self, tmp_path):
        content = "chrA\t500\nchrB\t700\n"
        fpath = tmp_path / "test.lens"
        fpath.write_text(content)
        data = read_lens(str(fpath))
        assert len(data) == 2

    def test_lens_all_filtered(self, tmp_path):
        lens_content = "chr1\t1000\nchr2\t2000\n"
        fpath = tmp_path / "test.lens"
        fpath.write_text(lens_content)

        chrs_content = "chrX\n"
        chrs_path = tmp_path / "chrs.txt"
        chrs_path.write_text(chrs_content)

        data = read_lens(str(fpath), chrs=str(chrs_path))
        assert len(data) == 0

    def test_single_chromosome(self, tmp_path):
        content = "chr1\t12345\n"
        fpath = tmp_path / "test.lens"
        fpath.write_text(content)
        data = read_lens(str(fpath))
        assert len(data) == 1
        assert data[0] == ["chr1", "12345"]


class TestProjectionInference:
    def test_brap_global_offset_prevents_local_chr4_flip(self):
        donors = [f"Chr{i}" for i in range(1, 5)]
        targets = [f"Chr{i}" for i in range(1, 9)]
        counts = {}
        for idx in range(1, 5):
            counts[(f"Chr{idx}", f"Chr{idx}")] = 20 if idx < 4 else 10
            counts[(f"Chr{idx}", f"Chr{idx + 4}")] = 10 if idx < 4 else 11
        raw, gff_a, gff_b, lens_a, lens_b = _projection_fixture(
            donors, targets, counts
        )

        interpreted, projection = infer_projection_color_map(
            raw, gff_a, gff_b, lens_a, lens_b, min_shared_pairs=5
        )

        assert projection[("Chr4", "Chr4")] == "primary"
        assert projection[("Chr4", "Chr8")] == "homeologous"
        assert interpreted["a_Chr4_1:b_Chr4_1"] == "red"
        assert interpreted["a_Chr4_1:b_Chr8_1"] == "blue"

    def test_bolr_global_second_offset_wins(self):
        donors = [f"Chr{i}" for i in range(1, 5)]
        targets = [f"Chr{i}" for i in range(1, 9)]
        counts = {}
        for idx in range(1, 5):
            counts[(f"Chr{idx}", f"Chr{idx}")] = 8
            counts[(f"Chr{idx}", f"Chr{idx + 4}")] = 20
        raw, gff_a, gff_b, lens_a, lens_b = _projection_fixture(
            donors, targets, counts
        )

        _, projection = infer_projection_color_map(
            raw, gff_a, gff_b, lens_a, lens_b, min_shared_pairs=5
        )

        assert projection[("Chr2", "Chr2")] == "homeologous"
        assert projection[("Chr2", "Chr6")] == "primary"

    @pytest.mark.parametrize(
        ("donor_suffix", "expected_primary"),
        [("A", "A"), ("D", "D")],
    )
    def test_wheat_diploid_suffix_mapping(self, donor_suffix, expected_primary):
        donors = [f"{idx}{donor_suffix}" for idx in range(1, 4)]
        targets = [f"{idx}{suffix}" for idx in range(1, 4) for suffix in "ABD"]
        counts = {}
        for idx in range(1, 4):
            for suffix in "ABD":
                counts[(f"{idx}{donor_suffix}", f"{idx}{suffix}")] = (
                    20 if suffix == donor_suffix else 8
                )
        raw, gff_a, gff_b, lens_a, lens_b = _projection_fixture(
            donors, targets, counts
        )

        _, projection = infer_projection_color_map(
            raw, gff_a, gff_b, lens_a, lens_b, min_shared_pairs=5
        )

        assert projection[(f"2{donor_suffix}", f"2{expected_primary}")] == "primary"
        other_suffix = "B" if donor_suffix != "B" else "D"
        assert projection[(f"2{donor_suffix}", f"2{other_suffix}")] == "homeologous"

    def test_wheat_ttur_preserves_a_and_b_suffixes(self):
        donors = [f"{idx}{suffix}" for idx in range(1, 3) for suffix in "AB"]
        targets = [f"{idx}{suffix}" for idx in range(1, 3) for suffix in "ABD"]
        counts = {}
        for idx in range(1, 3):
            for donor_suffix in "AB":
                for target_suffix in "ABD":
                    counts[(f"{idx}{donor_suffix}", f"{idx}{target_suffix}")] = (
                        20 if donor_suffix == target_suffix else 7
                    )
        raw, gff_a, gff_b, lens_a, lens_b = _projection_fixture(
            donors, targets, counts
        )

        _, projection = infer_projection_color_map(
            raw, gff_a, gff_b, lens_a, lens_b, min_shared_pairs=5
        )

        assert projection[("1A", "1A")] == "primary"
        assert projection[("1B", "1B")] == "primary"
        assert projection[("1A", "1D")] == "homeologous"
        assert projection[("1B", "1D")] == "homeologous"

    def test_close_offset_scores_are_ambiguous_and_keep_raw_color(self):
        donors = ["Chr1", "Chr2"]
        targets = ["Chr1", "Chr2", "Chr3", "Chr4"]
        counts = {
            ("Chr1", "Chr1"): 10, ("Chr2", "Chr2"): 10,
            ("Chr1", "Chr3"): 10, ("Chr2", "Chr4"): 10,
        }
        raw, gff_a, gff_b, lens_a, lens_b = _projection_fixture(
            donors, targets, counts
        )

        interpreted, projection = infer_projection_color_map(
            raw, gff_a, gff_b, lens_a, lens_b, min_shared_pairs=5
        )

        assert set(projection.values()) == {"ambiguous"}
        pair = "a_Chr1_1:b_Chr1_1"
        assert interpreted[pair] == raw[pair]
        assert confident_projection_pairs(projection) == []

    def test_one_to_many_pairs_do_not_overturn_unique_anchor_support(self):
        donors = ["Chr1", "Chr2"]
        targets = ["Chr1", "Chr2", "Chr3", "Chr4"]
        counts = {
            ("Chr1", "Chr1"): 20, ("Chr2", "Chr2"): 20,
            ("Chr1", "Chr3"): 8, ("Chr2", "Chr4"): 8,
        }
        raw, gff_a, gff_b, lens_a, lens_b = _projection_fixture(
            donors, targets, counts
        )
        for idx in range(20, 80):
            raw[f"a_Chr1_0:b_Chr3_{idx}"] = "blue"

        _, projection = infer_projection_color_map(
            raw, gff_a, gff_b, lens_a, lens_b, min_shared_pairs=5
        )

        assert projection[("Chr1", "Chr1")] == "primary"
        assert projection[("Chr1", "Chr3")] == "homeologous"

    def test_unstructured_and_missing_chromosomes_do_not_crash(self):
        raw, gff_a, gff_b, lens_a, lens_b = _projection_fixture(
            ["linkage_alpha", "linkage_beta"],
            ["group_x", "group_y", "group_unused"],
            {
                ("linkage_alpha", "group_x"): 20,
                ("linkage_beta", "group_y"): 20,
            },
        )
        raw["missing_donor:missing_target"] = "red"

        interpreted, projection = infer_projection_color_map(
            raw, gff_a, gff_b, lens_a, lens_b, min_shared_pairs=5
        )

        assert projection[("linkage_alpha", "group_x")] == "primary"
        assert projection[("linkage_beta", "group_y")] == "primary"
        assert interpreted["missing_donor:missing_target"] == "red"


class TestConversionZoneIndexing:
    def test_interleaved_chromosome_pairs_do_not_cross_contaminate_zones(self):
        data = []
        gff_a = {}
        gff_b = {}
        colors = ["red", "blue", "blue", "blue", "red"]
        for idx, color in enumerate(colors):
            for target_chr in ("Chr1", "Chr2"):
                gene_a = f"a_{target_chr}_{idx}"
                gene_b = f"b_{target_chr}_{idx}"
                data.append((gene_a, gene_b, color))
                gff_a[gene_a] = ["Chr1", str(idx + 1), str(idx + 1), "+", "."]
                gff_b[gene_b] = [target_chr, str(idx + 1), str(idx + 1), "+", "."]

        records = find_conversion_zones_with_ids_to_file(
            data, gff_a, gff_b, mode="intra", chrs_per_subgenome=10
        )

        by_zone = {}
        for record in records:
            by_zone.setdefault(record["zone_id"], []).append(record)
        assert len(by_zone) == 2
        assert all(len(zone_records) == 5 for zone_records in by_zone.values())
        for zone_records in by_zone.values():
            chromosome_pairs = {
                (record["contig1"], record["contig2"])
                for record in zone_records
            }
            assert len(chromosome_pairs) == 1


class TestSequenceSplitting:
    def test_complete_in_memory_assignments_are_not_limited_to_pair_tsv(self, tmp_path):
        fasta = tmp_path / "hybrid.fasta"
        fasta.write_text(">gene_a\nAAAA\n>gene_b\nCCCC\n>gene_c\nGGGG\n")
        gff = tmp_path / "hybrid.gff"
        gff.write_text(
            "Chr1\tgene_a\t1\t4\t+\t.\n"
            "Chr2\tgene_b\t1\t4\t+\t.\n"
            "Chr3\tgene_c\t1\t4\t+\t.\n"
        )

        split_sequences(
            str(fasta), str(tmp_path / "not_needed.tsv"), str(gff), str(tmp_path),
            subgenome_assignments={"gene_a": "A", "gene_b": "B", "gene_c": "unknown"},
        )

        summary = (tmp_path / "split_summary.txt").read_text()
        assert "assigned_total=2" in summary
        assert "skipped_unassigned=1" in summary
        assert ">gene_a" in (tmp_path / "split_subgenome_A.fasta").read_text()
        assert ">gene_b" in (tmp_path / "split_subgenome_B.fasta").read_text()
