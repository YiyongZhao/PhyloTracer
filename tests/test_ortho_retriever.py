"""
Tests for phylotracer.Ortho_Retriever module.
"""

import io
from itertools import combinations
import os
import sys
import textwrap
import tempfile
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    from ete3 import Tree
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False

from phylotracer.Ortho_Retriever import (
    rename_len_dic,
    count_sps_num,
    rm_dup,
    split_offcut_ev_seqs,
    parse_synteny_blocks,
    SyntenySource,
    filter_trees_by_synteny,
    choose_strict_single_copy_outgroup_genes,
    split_main,
)

EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), "..", "example_data", "14_Ortho_Retriever")
EXAMPLE_DIR = os.path.abspath(EXAMPLE_DIR)


# =============================================
# rename_len_dic
# =============================================

class TestRenameLenDic:
    def test_basic_rename(self):
        len_dic = {"geneA": 100, "geneB": 200, "geneC": 300}
        mapping = {"geneA": "sps_A", "geneB": "sps_B"}
        result = rename_len_dic(len_dic, mapping)
        assert result == {"sps_A": 100, "sps_B": 200}

    def test_no_overlap(self):
        result = rename_len_dic({"geneX": 50}, {"geneY": "sps_Y"})
        assert result == {}

    def test_empty_inputs(self):
        assert rename_len_dic({}, {}) == {}

    def test_full_overlap(self):
        len_dic = {"g1": 10, "g2": 20}
        mapping = {"g1": "r1", "g2": "r2"}
        result = rename_len_dic(len_dic, mapping)
        assert result == {"r1": 10, "r2": 20}


# =============================================
# count_sps_num
# =============================================

class TestCountSpsNum:
    def test_single_species(self):
        assert count_sps_num({"AAA_001", "AAA_002", "AAA_003"}) == 1

    def test_two_species(self):
        assert count_sps_num({"AAA_001", "BBB_001"}) == 2

    def test_empty(self):
        assert count_sps_num(set()) == 0

    def test_three_species(self):
        genes = {"AAA_001", "BBB_001", "CCC_001", "AAA_002"}
        assert count_sps_num(genes) == 3


# =============================================
# rm_dup
# =============================================

class TestRmDup:
    def test_no_subsets(self):
        a = {"g1", "g2"}
        b = {"g3", "g4"}
        result = rm_dup([a, b])
        assert len(result) == 2

    def test_proper_subset_removed(self):
        large = {"g1", "g2", "g3"}
        small = {"g1", "g2"}
        result = rm_dup([large, small])
        assert small not in result
        assert large in result

    def test_equal_sets_both_kept(self):
        a = {"g1", "g2"}
        b = {"g1", "g2"}
        result = rm_dup([a, b])
        assert len(result) == 2

    def test_empty_list(self):
        assert rm_dup([]) == []


# =============================================
# split_offcut_ev_seqs
# =============================================

class TestSplitOffcutEvSeqs:
    def test_ortholog_set(self):
        # 2 species, 2 genes — one gene per species → ortholog
        ortho = {"AAA_001", "BBB_001"}
        othologs_L, paralogs_L = split_offcut_ev_seqs([ortho])
        assert ortho in othologs_L
        assert paralogs_L == []

    def test_paralog_set(self):
        # 2 species, 3 genes → paralog
        para = {"AAA_001", "AAA_002", "BBB_001"}
        othologs_L, paralogs_L = split_offcut_ev_seqs([para])
        assert para in paralogs_L
        assert othologs_L == []

    def test_single_species_excluded(self):
        # only 1 species — should not appear in either list
        single = {"AAA_001", "AAA_002"}
        othologs_L, paralogs_L = split_offcut_ev_seqs([single])
        assert othologs_L == []
        assert paralogs_L == []

    def test_mixed(self):
        ortho = {"AAA_001", "BBB_001"}
        para = {"AAA_001", "AAA_002", "BBB_001"}
        othologs_L, paralogs_L = split_offcut_ev_seqs([ortho, para])
        assert ortho in othologs_L
        assert para in paralogs_L


# =============================================
# parse_synteny_blocks
# =============================================

class TestParseSyntenyBlocks:
    def _write_tmp(self, content):
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False, encoding="utf-8")
        f.write(textwrap.dedent(content))
        f.flush()
        f.close()
        return f.name

    def test_basic_two_column(self):
        path = self._write_tmp("""\
            # block header
            geneA\tgeneB
            geneC\tgeneD
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert result["geneA"] == {"geneB": 1.0}
        assert result["geneB"] == {"geneA": 1.0}

    def test_wgdi_style(self):
        # gene1  order  gene2  order  strand
        path = self._write_tmp("""\
            # Alignment\t1:
            geneA\t1\tgeneB\t2\t1
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert result["geneA"] == {"geneB": 1.0}
        assert result["geneB"] == {"geneA": 1.0}

    def test_empty_file(self):
        path = self._write_tmp("")
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert result == {}

    def test_multiple_blocks(self):
        path = self._write_tmp("""\
            # block1
            g1\tg2
            # block2
            g3\tg4
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert result["g1"] == {"g2": 1.0}
        assert result["g3"] == {"g4": 1.0}

    def test_gene_in_multiple_blocks(self):
        path = self._write_tmp("""\
            # block1
            shared\tg2
            # block2
            shared\tg3
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert result["shared"] == {"g2": 1.0, "g3": 1.0}

    def test_three_column_anchors_and_reversed_duplicates(self):
        path = self._write_tmp("""\
            ###
            AAB_1\tAAC_1\t131
            AAC_1\tAAB_1\t2140
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert result == {"AAB_1": {"AAC_1": 2140.0}, "AAC_1": {"AAB_1": 2140.0}}

    def test_duplicate_anchor_score_uses_last_row(self):
        path = self._write_tmp("""\
            AAB_1\tAAC_1\t2140
            AAC_1\tAAB_1\t131
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert result["AAB_1"]["AAC_1"] == 131.0

    def test_genes_in_same_block_are_not_implicitly_paired(self):
        path = self._write_tmp("""\
            # block1
            A_1\tB_1
            A_2\tB_2
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert "B_2" not in result["A_1"]


# =============================================
# filter_trees_by_synteny
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestFilterTreesBySynteny:
    def _make_tree(self, leaf_names):
        newick = "(" + ",".join(f"{n}:1" for n in leaf_names) + ");"
        return Tree(newick)

    def test_missing_pairs_count_in_denominator(self):
        tree = self._make_tree(["g1", "g2", "g3"])
        species = {"g1": "A", "g2": "B", "g3": "C"}
        filtered, report = filter_trees_by_synteny(
            [("T1", tree)], {"g1": {"g2"}}, species, synteny_threshold=0.3
        )
        assert len(filtered) == 1
        assert report[0][2:7] == ("3", "1", "2", "0", "0.333333")
        assert report[0][-1] == "pass"

    def test_fails_with_no_support(self):
        tree = self._make_tree(["g1", "g2", "g3"])
        species = {"g1": "A", "g2": "B", "g3": "C"}
        filtered, report = filter_trees_by_synteny([("T1", tree)], {}, species)
        assert len(filtered) == 0
        assert report[0][2:7] == ("3", "0", "3", "0", "0.000000")
        assert report[0][-1] == "fail"

    def test_threshold_boundary_uses_exact_rate(self):
        tree = self._make_tree(["g1", "g2", "g3", "g4"])
        species = {gene: gene for gene in ["g1", "g2", "g3", "g4"]}
        partners = {"g1": {"g2", "g3"}}
        passed, report = filter_trees_by_synteny(
            [("T1", tree)], partners, species, synteny_threshold=1 / 3
        )
        failed, _ = filter_trees_by_synteny(
            [("T1", tree)], partners, species, synteny_threshold=0.34
        )
        assert len(passed) == 1
        assert len(failed) == 0
        assert report[0][2:5] == ("6", "2", "4")

    def test_duplicate_species_rejected_as_not_single_copy(self):
        tree = self._make_tree(["g1", "g2", "g3"])
        species = {"g1": "A", "g2": "A", "g3": "B"}
        partners = {"g1": {"g2", "g3"}}
        with pytest.raises(ValueError, match="not single-copy"):
            filter_trees_by_synteny(
                [("T1", tree)], partners, species, synteny_threshold=0.5
            )

    def test_one_tip_fails_even_at_zero_threshold(self):
        tree = Tree("g1;")
        filtered, report = filter_trees_by_synteny(
            [("T1", tree)], {}, {"g1": "A"},
            synteny_threshold=0,
        )
        assert report[0][-1] == "fail"
        assert filtered == []

    def test_missing_species_mapping_fails_explicitly(self):
        tree = self._make_tree(["g1", "g2"])
        with pytest.raises(ValueError, match="Missing species mapping"):
            filter_trees_by_synteny([("T1", tree)], {}, {"g1": "A"})

    def test_report_columns(self):
        tree = self._make_tree(["g1", "g2"])
        _, report = filter_trees_by_synteny(
            [("T1", tree)], {"g1": {"g2"}}, {"g1": "A", "g2": "B"}
        )
        row = report[0]
        assert len(row) == 10
        assert row[0] == "T1"

    def test_missing_species_file_stays_in_denominator(self, tmp_path):
        (tmp_path / "A.B.anchors").write_text("###\nA_1\tB_1\t131\n")
        (tmp_path / "A.C.anchors").write_text("###\nA_2\tC_2\t2140\n")
        source = SyntenySource(str(tmp_path), {"A_1": "A", "B_1": "B", "C_1": "C"})
        tree = self._make_tree(["A_1", "B_1", "C_1"])
        filtered, report = filter_trees_by_synteny(
            [("T1", tree)], source,
            {"A_1": "A", "B_1": "B", "C_1": "C"},
            synteny_threshold=0.3,
        )
        assert len(filtered) == 1
        assert report[0][2:8] == ("3", "1", "2", "1", "0.333333", "131.00")
        failed, _ = filter_trees_by_synteny(
            [("T1", tree)], source,
            {"A_1": "A", "B_1": "B", "C_1": "C"},
            synteny_threshold=0.4,
        )
        assert failed == []

    def test_empty_species_pair_file_is_not_assessable(self, tmp_path):
        (tmp_path / "A.B.anchors").write_text("###\n")
        source = SyntenySource(str(tmp_path), {"A_1": "A", "B_1": "B"})
        tree = self._make_tree(["A_1", "B_1"])
        filtered, report = filter_trees_by_synteny(
            [("T1", tree)], source, {"A_1": "A", "B_1": "B"},
            synteny_threshold=0.2,
        )
        assert filtered == []
        assert report[0][2:7] == ("1", "0", "1", "1", "0.000000")

    def test_directory_index_reuses_and_refreshes_anchor_pairs(self, tmp_path, monkeypatch):
        monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
        anchor_dir = tmp_path / "anchors"
        anchor_dir.mkdir()
        anchor = anchor_dir / "A.B.anchors"
        anchor.write_text("A_1\tB_1\t131\nB_1\tA_1\t2140\n")
        species = {"A_1": "A", "B_1": "B"}

        source = SyntenySource(str(anchor_dir), species)
        assert source.pair_score("A_1", "B_1", "A", "B") == (True, 2140.0)
        assert source.pair_score("B_1", "A_1", "B", "A") == (True, 2140.0)
        assert source.pair_score("A_2", "B_2", "A", "B") == (True, None)
        source.close()

        from phylotracer import Ortho_Retriever as module
        original_iterator = module._iter_synteny_pairs
        def fail_if_reparsed(_):
            raise AssertionError("unchanged anchors file was reparsed")
        monkeypatch.setattr(module, "_iter_synteny_pairs", fail_if_reparsed)
        source = SyntenySource(str(anchor_dir), species)
        assert source.pair_score("A_1", "B_1", "A", "B") == (True, 2140.0)
        source.close()
        monkeypatch.setattr(module, "_iter_synteny_pairs", original_iterator)

        anchor.write_text("A_1\tB_1\t42\n")
        source = SyntenySource(str(anchor_dir), species)
        assert source.pair_score("A_1", "B_1", "A", "B") == (True, 42.0)
        source.close()

    def test_ten_tips_twenty_seven_matches_use_all_forty_five_pairs(self):
        genes = [f"S{i}_g" for i in range(10)]
        species = {gene: gene.split("_", 1)[0] for gene in genes}
        partners = {}
        for gene1, gene2 in list(combinations(genes, 2))[:27]:
            partners.setdefault(gene1, set()).add(gene2)
        tree = self._make_tree(genes)

        passed, report = filter_trees_by_synteny(
            [("T10", tree)], partners, species, synteny_threshold=0.6
        )
        failed, _ = filter_trees_by_synteny(
            [("T10", tree)], partners, species, synteny_threshold=0.61
        )

        assert len(passed) == 1
        assert failed == []
        assert report[0][2:7] == ("45", "27", "18", "0", "0.600000")

    @pytest.mark.parametrize("threshold", [-0.01, 1.01])
    def test_invalid_threshold_rejected(self, threshold):
        with pytest.raises(ValueError, match="synteny_threshold"):
            filter_trees_by_synteny([], {}, {}, synteny_threshold=threshold)


# =============================================
# strict outgroup helpers
# =============================================

class TestStrictOutgroupSelection:
    def test_filters_overlap_and_multicopy(self):
        selected = choose_strict_single_copy_outgroup_genes(
            candidate_genes=["g1", "g2", "g3"],
            gene_to_species={"g1": "sp1", "g2": "sp2", "g3": "sp2"},
            ingroup_species={"sp1"},
        )
        assert selected == []

    def test_keeps_single_copy_nonoverlap(self):
        selected = choose_strict_single_copy_outgroup_genes(
            candidate_genes=["g1", "g2", "g3"],
            gene_to_species={"g1": "sp1", "g2": "sp2", "g3": "sp3"},
            ingroup_species={"sp1"},
        )
        assert selected == ["g2", "g3"]


# =============================================
# Integration: split_main with example data
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestSplitMainIntegration:
    """Run split_main on the bundled example_data and verify outputs."""

    @pytest.fixture(autouse=True)
    def run_in_tmp(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        self.tmp = tmp_path

    def _load_imap(self, path):
        d = {}
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) == 2:
                    d[parts[0]] = parts[1]
        return d

    def test_output_files_created(self):
        from phylotracer import gene_id_transfer, read_and_return_dict
        imap_path = os.path.join(EXAMPLE_DIR, "gene2sps.imap")
        len_path = os.path.join(EXAMPLE_DIR, "gene2length.imap")
        gf_path = os.path.join(EXAMPLE_DIR, "GF_ID2path.imap")

        gene2new, new2gene, _, _ = gene_id_transfer(imap_path)
        len_dic = read_and_return_dict(len_path)
        renamed_len_dic = rename_len_dic(len_dic, gene2new)

        # GF_ID2path.imap uses relative paths from example_data directory
        raw_tre_dic = read_and_return_dict(gf_path)
        base = os.path.join(EXAMPLE_DIR)
        tre_dic = {
            k: os.path.normpath(os.path.join(base, v))
            for k, v in raw_tre_dic.items()
        }

        split_main(tre_dic, gene2new, new2gene, renamed_len_dic)

        summary = self.tmp / "ortho_retriever_summary.txt"
        tsv = self.tmp / "ortholog_trees.tsv"
        assert summary.exists(), "ortho_retriever_summary.txt not created"
        assert tsv.exists(), "ortholog_trees.tsv not created"

    def test_summary_has_header_and_data(self):
        from phylotracer import gene_id_transfer, read_and_return_dict
        imap_path = os.path.join(EXAMPLE_DIR, "gene2sps.imap")
        len_path = os.path.join(EXAMPLE_DIR, "gene2length.imap")
        gf_path = os.path.join(EXAMPLE_DIR, "GF_ID2path.imap")

        gene2new, new2gene, _, _ = gene_id_transfer(imap_path)
        len_dic = read_and_return_dict(len_path)
        renamed_len_dic = rename_len_dic(len_dic, gene2new)
        raw_tre_dic = read_and_return_dict(gf_path)
        base = os.path.join(EXAMPLE_DIR)
        tre_dic = {
            k: os.path.normpath(os.path.join(base, v))
            for k, v in raw_tre_dic.items()
        }

        split_main(tre_dic, gene2new, new2gene, renamed_len_dic)

        lines = (self.tmp / "ortho_retriever_summary.txt").read_text().splitlines()
        assert lines[0] == "tre_name\tsingle_copy_tree"
        assert len(lines) > 1, "summary file has no data rows"

    def test_synteny_filter_produces_report(self):
        from phylotracer import gene_id_transfer, read_and_return_dict
        imap_path = os.path.join(EXAMPLE_DIR, "gene2sps.imap")
        len_path = os.path.join(EXAMPLE_DIR, "gene2length.imap")
        gf_path = os.path.join(EXAMPLE_DIR, "GF_ID2path.imap")
        synteny_path = os.path.join(EXAMPLE_DIR, "collinearity")

        gene2new, new2gene, _, _ = gene_id_transfer(imap_path)
        len_dic = read_and_return_dict(len_path)
        renamed_len_dic = rename_len_dic(len_dic, gene2new)
        raw_tre_dic = read_and_return_dict(gf_path)
        base = os.path.join(EXAMPLE_DIR)
        tre_dic = {
            k: os.path.normpath(os.path.join(base, v))
            for k, v in raw_tre_dic.items()
        }

        split_main(
            tre_dic, gene2new, new2gene, renamed_len_dic,
            synteny_blocks_path=synteny_path,
            gene_to_species=read_and_return_dict(imap_path),
        )

        report = self.tmp / "ortholog_synteny_report.tsv"
        assert report.exists(), "ortholog_synteny_report.tsv not created"
        lines = report.read_text().splitlines()
        assert lines[0] == (
            "tre_ID\ttre_name\tnum_genes\ttotal_pairs\tmatched_pairs\t"
            "unmatched_pairs\tpairs_without_synteny_data\tsynteny_rate\t"
            "mean_matched_anchor_score\tthreshold\tsynteny_filter"
        )
        assert len(lines) > 1

    def test_all_filtered_out_replaces_stale_ortholog_table(self):
        from phylotracer import gene_id_transfer, read_and_return_dict

        imap_path = os.path.join(EXAMPLE_DIR, "gene2sps.imap")
        len_path = os.path.join(EXAMPLE_DIR, "gene2length.imap")
        gf_path = os.path.join(EXAMPLE_DIR, "GF_ID2path.imap")
        gene2new, new2gene, _, _ = gene_id_transfer(imap_path)
        len_dic = read_and_return_dict(len_path)
        raw_tre_dic = read_and_return_dict(gf_path)
        first_id, first_path = next(iter(raw_tre_dic.items()))
        tre_dic = {first_id: os.path.join(EXAMPLE_DIR, first_path)}
        anchors = self.tmp / "empty.anchors"
        anchors.write_text("###\n")
        (self.tmp / "ortholog_trees.tsv").write_text("stale\n")
        (self.tmp / "ortholog_outgroup_report.tsv").write_text("stale\n")

        split_main(
            tre_dic, gene2new, new2gene, rename_len_dic(len_dic, gene2new),
            synteny_blocks_path=str(anchors),
            gene_to_species=read_and_return_dict(imap_path),
            add_outgroup=True,
        )

        assert (self.tmp / "ortholog_trees.tsv").read_text().strip() == "tre_ID"
        assert (self.tmp / "ortholog_outgroup_report.tsv").read_text().splitlines() == [
            "tre_name\tgf_id\tnum_ingroup\tnum_outgroup\toutgroup_genes\tstatus"
        ]

    def test_add_outgroup_produces_report(self):
        from phylotracer import gene_id_transfer, read_and_return_dict

        imap_path = os.path.join(EXAMPLE_DIR, "gene2sps.imap")
        len_path = os.path.join(EXAMPLE_DIR, "gene2length.imap")
        gf_path = os.path.join(EXAMPLE_DIR, "GF_ID2path.imap")

        gene_to_species = read_and_return_dict(imap_path)
        gene2new, new2gene, _, _ = gene_id_transfer(imap_path)
        len_dic = read_and_return_dict(len_path)
        renamed_len_dic = rename_len_dic(len_dic, gene2new)
        raw_tre_dic = read_and_return_dict(gf_path)
        base = os.path.join(EXAMPLE_DIR)
        tre_dic = {
            k: os.path.normpath(os.path.join(base, v))
            for k, v in raw_tre_dic.items()
        }

        split_main(
            tre_dic,
            gene2new,
            new2gene,
            renamed_len_dic,
            add_outgroup=True,
            gene_to_species=gene_to_species,
        )

        report = self.tmp / "ortholog_outgroup_report.tsv"
        assert report.exists(), "ortholog_outgroup_report.tsv not created"
        lines = report.read_text().splitlines()
        assert lines[0].startswith("tre_name\tgf_id\t")
        assert len(lines) > 1
        assert any("\tok" in line or "\tskip_no_valid_outgroup" in line for line in lines[1:])
