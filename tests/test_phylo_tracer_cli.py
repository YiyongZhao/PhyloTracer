"""
Tests for phylotracer.Phylo_Tracer CLI entry point.
"""

import os
import sys
import subprocess
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
            cwd="/tmp/PhyloTracer",
        )
        # argparse -h exits with code 0
        assert result.returncode == 0

    def test_help_output_contains_phylotracer(self):
        result = subprocess.run(
            [sys.executable, "-m", "phylotracer.Phylo_Tracer", "-h"],
            capture_output=True,
            text=True,
            cwd="/tmp/PhyloTracer",
        )
        assert "PhyloTracer" in result.stdout or "phylotracer" in result.stdout.lower()

    def test_help_lists_subcommands(self):
        result = subprocess.run(
            [sys.executable, "-m", "phylotracer.Phylo_Tracer", "-h"],
            capture_output=True,
            text=True,
            cwd="/tmp/PhyloTracer",
        )
        assert "GD_Detector" in result.stdout
        assert "Phylo_Rooter" in result.stdout


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
            cwd="/tmp/PhyloTracer",
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
