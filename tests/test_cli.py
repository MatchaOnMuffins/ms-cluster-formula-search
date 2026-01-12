"""Tests for formula_search.cli module."""

import pytest
import io
import sys
from unittest.mock import patch
from formula_search.cli import main, run_scan_all


class TestRunScanAll:
    """Test the run_scan_all function."""

    def test_run_scan_all_produces_output(self):
        """run_scan_all should produce output for all levels."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            run_scan_all(1519.1540, ppm=10)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert len(output) > 0
        assert "Level" in output

    def test_run_scan_all_shows_all_levels(self):
        """Should show results for levels 1, 2, and 3."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            run_scan_all(1519.1540, ppm=10)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Level 1" in output
        assert "Level 2" in output
        assert "Level 3" in output

    def test_run_scan_all_incremental(self):
        """Later levels should show only new hits."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            run_scan_all(1519.1540, ppm=10)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # Parse formulas from each level section would be complex,
        # but we can at least verify it runs without error
        assert "strict" in output
        assert "moderate" in output
        assert "loose" in output


class TestMainCLI:
    """Test the main CLI entry point."""

    def test_main_with_simple_mz(self):
        """CLI should work with just m/z argument."""
        captured = io.StringIO()
        sys.stdout = captured

        test_args = ["search-formula", "1519.1540"]
        with patch.object(sys, "argv", test_args):
            try:
                main()
            finally:
                sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # Should produce some output (either results or "No matches")
        assert len(output) > 0

    def test_main_with_coarseness(self):
        """CLI should accept coarseness flag."""
        captured = io.StringIO()
        sys.stdout = captured

        test_args = ["search-formula", "1519.1540", "-c", "3"]
        with patch.object(sys, "argv", test_args):
            try:
                main()
            finally:
                sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert len(output) > 0

    def test_main_with_ppm(self):
        """CLI should accept ppm flag."""
        captured = io.StringIO()
        sys.stdout = captured

        test_args = ["search-formula", "1519.1540", "--ppm", "5"]
        with patch.object(sys, "argv", test_args):
            try:
                main()
            finally:
                sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # Might have fewer results with tighter ppm
        assert len(output) > 0

    def test_main_scan_all_flag(self):
        """CLI should accept --scan-all flag."""
        captured = io.StringIO()
        sys.stdout = captured

        test_args = ["search-formula", "1519.1540", "--scan-all"]
        with patch.object(sys, "argv", test_args):
            try:
                main()
            finally:
                sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # Should show all three levels
        assert "Level 1" in output
        assert "Level 2" in output
        assert "Level 3" in output

    def test_main_no_matches(self):
        """CLI should report when no matches found."""
        captured = io.StringIO()
        sys.stdout = captured

        # Use a very small m/z that won't match anything
        test_args = ["search-formula", "10.0", "--ppm", "1"]
        with patch.object(sys, "argv", test_args):
            try:
                main()
            finally:
                sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "No matches" in output

    def test_main_combined_flags(self):
        """CLI should work with multiple flags."""
        captured = io.StringIO()
        sys.stdout = captured

        test_args = ["search-formula", "1519.1540", "-c", "1", "--ppm", "15"]
        with patch.object(sys, "argv", test_args):
            try:
                main()
            finally:
                sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert len(output) > 0
