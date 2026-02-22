"""
Tests for phylotracer.HaploFinder module.
"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from phylotracer.HaploFinder import read_gff, read_lens


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
