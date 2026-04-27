"""
Tests for phylotracer.Ortho_Retriever module.
"""

import io
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
        assert "geneA" in result
        assert "geneB" in result

    def test_wgdi_style(self):
        # gene1  order  gene2  order  strand
        path = self._write_tmp("""\
            # Alignment\t1:
            geneA\t1\tgeneB\t2\t1
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert "geneA" in result
        assert "geneB" in result

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
        assert len(result["g1"]) == 1
        assert len(result["g3"]) == 1
        # g1 and g3 should belong to different blocks
        assert result["g1"] != result["g3"]

    def test_gene_in_multiple_blocks(self):
        path = self._write_tmp("""\
            # block1
            shared\tg2
            # block2
            shared\tg3
        """)
        result = parse_synteny_blocks(path)
        os.unlink(path)
        assert len(result["shared"]) == 2


# =============================================
# filter_trees_by_synteny
# =============================================

@pytest.mark.skipif(not HAS_ETE3, reason="ete3 not installed")
class TestFilterTreesBySynteny:
    def _make_tree(self, leaf_names):
        newick = "(" + ",".join(f"{n}:1" for n in leaf_names) + ");"
        return Tree(newick)

    def test_passes_with_enough_support(self):
        tree = self._make_tree(["g1", "g2", "g3"])
        gene_to_blocks = {
            "g1": {"blk1"},
            "g2": {"blk1"},
            "g3": {"blk1"},
        }
        filtered, report = filter_trees_by_synteny([("T1", tree)], gene_to_blocks)
        assert len(filtered) == 1
        assert report[0][-1] == "pass"

    def test_fails_with_no_support(self):
        tree = self._make_tree(["g1", "g2", "g3"])
        filtered, report = filter_trees_by_synteny([("T1", tree)], {})
        assert len(filtered) == 0
        assert report[0][-1] == "fail"

    def test_fails_when_support_ratio_low(self):
        # g1 in blk1, g2 in blk2 → best_block covers 1/2 genes, ratio=0.5 → pass threshold exactly
        tree = self._make_tree(["g1", "g2"])
        gene_to_blocks = {"g1": {"blk1"}, "g2": {"blk2"}}
        filtered, report = filter_trees_by_synteny([("T1", tree)], gene_to_blocks)
        # supported_genes=2, best_block_gene_count=1 → ratio=0.5 → passes (>=0.5)
        assert report[0][-1] == "pass"

    def test_fails_when_only_one_supported_gene(self):
        tree = self._make_tree(["g1", "g2", "g3"])
        gene_to_blocks = {"g1": {"blk1"}}
        filtered, report = filter_trees_by_synteny([("T1", tree)], gene_to_blocks)
        assert report[0][-1] == "fail"

    def test_report_columns(self):
        tree = self._make_tree(["g1", "g2"])
        gene_to_blocks = {"g1": {"blk1"}, "g2": {"blk1"}}
        _, report = filter_trees_by_synteny([("T1", tree)], gene_to_blocks)
        row = report[0]
        # (tree_name, num_genes, supported_genes, best_block, best_block_gene_count, ratio, verdict)
        assert len(row) == 7
        assert row[0] == "T1"


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

        split_main(tre_dic, gene2new, new2gene, renamed_len_dic, synteny_blocks_path=synteny_path)

        report = self.tmp / "ortholog_synteny_report.tsv"
        assert report.exists(), "ortholog_synteny_report.tsv not created"
        lines = report.read_text().splitlines()
        assert lines[0].startswith("tre_ID\t")
        assert len(lines) > 1

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
