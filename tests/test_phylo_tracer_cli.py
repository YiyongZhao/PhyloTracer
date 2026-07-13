"""
Tests for phylotracer.Phylo_Tracer CLI entry point.
"""

import os
import sys
import subprocess
from types import SimpleNamespace
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


# =============================================
# main function existence tests
# =============================================

class TestMainExists:
    def test_main_is_callable(self):
        from phylotracer.Phylo_Tracer import main
        assert callable(main)

    def test_command_handlers_is_dict(self):
        from phylotracer.Phylo_Tracer import COMMAND_HANDLERS
        assert isinstance(COMMAND_HANDLERS, dict)

    def test_command_handlers_has_entries(self):
        from phylotracer.Phylo_Tracer import COMMAND_HANDLERS
        assert len(COMMAND_HANDLERS) > 0

    def test_all_expected_subcommands_present(self):
        from phylotracer.Phylo_Tracer import COMMAND_HANDLERS
        expected = [
            "PhyloTree_CollapseExpand",
            "Phylo_Rooter",
            "OrthoFilter_LB",
            "OrthoFilter_Mono",
            "TreeTopology_Summarizer",
            "GD_Detector",
            "HaploFinder",
        ]
        for cmd in expected:
            assert cmd in COMMAND_HANDLERS, f"Missing subcommand: {cmd}"

    def test_parser_is_available(self):
        from phylotracer.Phylo_Tracer import _parser
        assert _parser is not None


# =============================================
# CLI help output tests
# =============================================

class TestCLIHelp:
    def test_help_flag_exits_zero(self):
        result = subprocess.run(
            [sys.executable, "-m", "phylotracer.Phylo_Tracer", "-h"],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        )
        # argparse -h exits with code 0
        assert result.returncode == 0

    def test_help_output_contains_phylotracer(self):
        result = subprocess.run(
            [sys.executable, "-m", "phylotracer.Phylo_Tracer", "-h"],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        )
        assert "PhyloTracer" in result.stdout or "phylotracer" in result.stdout.lower()

    def test_help_lists_subcommands(self):
        result = subprocess.run(
            [sys.executable, "-m", "phylotracer.Phylo_Tracer", "-h"],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        )
        assert "GD_Detector" in result.stdout
        assert "Phylo_Rooter" in result.stdout

    def test_haplofinder_help_hides_legacy_no_effect_options(self):
        result = subprocess.run(
            [sys.executable, "-m", "phylotracer.Phylo_Tracer", "HaploFinder", "--help"],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        )
        assert result.returncode == 0
        assert "--n_permutations" not in result.stdout
        assert "--p_threshold" not in result.stdout
        assert "--cluster_file" not in result.stdout


class TestHaploFinderCLI:
    required_args = [
        "--input_GF_list", "gf.imap",
        "--input_imap", "genes.imap",
        "--input_sps_tree", "species.nwk",
        "--species_a", "ARD",
        "--species_b", "ARH",
        "--species_a_gff", "ARD.gff",
        "--species_b_gff", "ARH.gff",
        "--species_a_lens", "ARD.lens",
        "--species_b_lens", "ARH.lens",
    ]

    def test_missing_core_arguments_exits_two(self):
        from phylotracer.Phylo_Tracer import _parser
        with pytest.raises(SystemExit) as exc_info:
            _parser.parse_args(["HaploFinder"])
        assert exc_info.value.code == 2

    def test_default_command_needs_no_projection_parameter(self):
        from phylotracer.Phylo_Tracer import _parser
        args = _parser.parse_args(["HaploFinder", *self.required_args])
        assert args.min_shared_pairs == 5
        assert not hasattr(args, "projection_mode")

    def test_legacy_arguments_remain_parseable(self):
        from phylotracer.Phylo_Tracer import _parser
        args = _parser.parse_args([
            "HaploFinder", *self.required_args,
            "--n_permutations", "200",
            "--p_threshold", "0.01",
            "--cluster_file", "legacy.tsv",
        ])
        assert args.n_permutations == 200
        assert args.p_threshold == pytest.approx(0.01)
        assert args.cluster_file == "legacy.tsv"

    @pytest.mark.parametrize(
        "extra_args",
        [
            ["--hyb_sps", "ARH"],
            ["--parental_sps", "ARD ARI"],
            ["--input_fasta", "ARH.pep"],
        ],
    )
    def test_subgenome_dependency_validation_exits_two(self, extra_args):
        result = subprocess.run(
            [
                sys.executable, "-m", "phylotracer.Phylo_Tracer", "HaploFinder",
                *self.required_args, *extra_args,
            ],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        )
        assert result.returncode == 2

    def test_legacy_arguments_emit_deprecation_warning(self, monkeypatch, caplog, tmp_path):
        import phylotracer.Phylo_Tracer as cli
        monkeypatch.setattr(cli, "process_gd_result", lambda *args, **kwargs: ([], {}, {}))
        monkeypatch.setattr(cli, "generate_dotplot", lambda *args, **kwargs: None)
        monkeypatch.setattr(cli, "report_execution_time", lambda *args, **kwargs: None)
        args = SimpleNamespace(
            input_GF_list="gf.imap", input_imap="genes.imap",
            input_sps_tree="species.nwk", species_a="ARD", species_b="ARH",
            species_a_gff="ARD.gff", species_b_gff="ARH.gff",
            species_a_lens="ARD.lens", species_b_lens="ARH.lens",
            gd_support=50, pair_support=50, size=0.0005,
            hyb_sps=None, parental_sps=None, input_fasta=None,
            visual_chr_a=None, visual_chr_b=None, min_conv_pairs=10,
            min_shared_pairs=5, output_dir=str(tmp_path),
            n_permutations=200, p_threshold=0.01, cluster_file="legacy.tsv",
        )

        with caplog.at_level("WARNING"):
            cli.handle_haplofinder(args)

        assert "--n_permutations is deprecated and has no effect" in caplog.text
        assert "--p_threshold is deprecated and has no effect" in caplog.text
        assert "--cluster_file is deprecated and has no effect" in caplog.text


# =============================================
# Invalid command handling tests
# =============================================

class TestInvalidCommand:
    def test_no_command_prints_usage(self):
        """Running without a subcommand should print usage and exit."""
        from phylotracer.Phylo_Tracer import main
        # main() with no args prints usage and returns None
        # We test indirectly: _parser.parse_args([]) should set command=None
        from phylotracer.Phylo_Tracer import _parser
        args = _parser.parse_args([])
        assert args.command is None

    def test_invalid_subcommand_causes_error(self):
        result = subprocess.run(
            [sys.executable, "-m", "phylotracer.Phylo_Tracer", "NonExistentCommand"],
            capture_output=True,
            text=True,
            cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        )
        # argparse should exit non-zero for invalid subcommand
        assert result.returncode != 0

    def test_print_cli_usage_callable(self):
        from phylotracer.Phylo_Tracer import print_cli_usage
        assert callable(print_cli_usage)


# =============================================
# Bounded validators tests
# =============================================

class TestBoundedValidators:
    def test_bounded_int_valid(self):
        from phylotracer.Phylo_Tracer import bounded_int
        validator = bounded_int(0, 100)
        assert validator("50") == 50

    def test_bounded_int_too_low(self):
        from phylotracer.Phylo_Tracer import bounded_int
        import argparse
        validator = bounded_int(0, 100)
        with pytest.raises(argparse.ArgumentTypeError):
            validator("-1")

    def test_bounded_int_too_high(self):
        from phylotracer.Phylo_Tracer import bounded_int
        import argparse
        validator = bounded_int(0, 100)
        with pytest.raises(argparse.ArgumentTypeError):
            validator("101")

    def test_bounded_float_valid(self):
        from phylotracer.Phylo_Tracer import bounded_float
        validator = bounded_float(0.0, 1.0)
        assert validator("0.5") == pytest.approx(0.5)

    def test_bounded_float_too_low(self):
        from phylotracer.Phylo_Tracer import bounded_float
        import argparse
        validator = bounded_float(0.0, 1.0)
        with pytest.raises(argparse.ArgumentTypeError):
            validator("-0.1")

    def test_bounded_float_too_high(self):
        from phylotracer.Phylo_Tracer import bounded_float
        import argparse
        validator = bounded_float(0.0, 1.0)
        with pytest.raises(argparse.ArgumentTypeError):
            validator("1.5")

    def test_bounded_int_boundary_values(self):
        from phylotracer.Phylo_Tracer import bounded_int
        validator = bounded_int(0, 100)
        assert validator("0") == 0
        assert validator("100") == 100

    def test_bounded_float_boundary_values(self):
        from phylotracer.Phylo_Tracer import bounded_float
        validator = bounded_float(0.0, 1.0)
        assert validator("0.0") == 0.0
        assert validator("1.0") == 1.0
