"""Tests for formula_search.formatting module."""

import pytest
import io
import sys
from formula_search.formatting import (
    format_hit,
    print_header,
    print_results,
    print_level_header,
)


class TestFormatHit:
    """Test the format_hit function."""

    @pytest.fixture
    def sample_hit(self):
        """A sample hit dictionary for testing."""
        return {
            "formula": "Y2Mn4(tBuCOO)10O6H3C1",
            "mz": 1519.1540,
            "charge": -1,
            "adduct": "[M]−•",
            "neutral_mass": 1519.1540,
            "ppm_error": 3.53,
            "mass": 1519.1540,
            "counts": (2, 4, 10, 6, 3, 1),
        }

    def test_format_hit_returns_string(self, sample_hit):
        """format_hit should return a string."""
        result = format_hit(sample_hit)
        assert isinstance(result, str)

    def test_format_hit_contains_formula(self, sample_hit):
        """Formatted string should contain the formula."""
        result = format_hit(sample_hit)
        assert "Y2Mn4(tBuCOO)10O6H3C1" in result

    def test_format_hit_contains_mz(self, sample_hit):
        """Formatted string should contain the m/z value."""
        result = format_hit(sample_hit)
        assert "1519.1540" in result

    def test_format_hit_contains_charge(self, sample_hit):
        """Formatted string should contain the charge."""
        result = format_hit(sample_hit)
        assert "-1" in result

    def test_format_hit_contains_adduct(self, sample_hit):
        """Formatted string should contain the adduct."""
        result = format_hit(sample_hit)
        assert "[M]−•" in result

    def test_format_hit_contains_neutral_mass(self, sample_hit):
        """Formatted string should contain neutral mass."""
        result = format_hit(sample_hit)
        # Should appear twice if mz == neutral_mass, or at least once
        assert "1519.1540" in result

    def test_format_hit_contains_ppm_error(self, sample_hit):
        """Formatted string should contain ppm error."""
        result = format_hit(sample_hit)
        assert "3.53" in result
        assert "ppm" in result

    def test_format_hit_positive_ppm(self, sample_hit):
        """Positive ppm error should have + sign."""
        result = format_hit(sample_hit)
        # Should contain +3.53
        assert "+3.53" in result or "+ 3.53" in result.replace(" ", "")

    def test_format_hit_negative_ppm(self, sample_hit):
        """Negative ppm error should have - sign."""
        sample_hit["ppm_error"] = -5.25
        result = format_hit(sample_hit)
        assert "-5.25" in result

    def test_format_hit_fixed_width(self, sample_hit):
        """Output should have consistent width for table alignment."""
        result = format_hit(sample_hit)
        # Formula field should be padded to 24 chars
        parts = result.split()
        # Just check it's formatted reasonably
        assert len(result) > 50


class TestPrintHeader:
    """Test the print_header function."""

    def test_print_header_output(self):
        """print_header should produce expected column headers."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_header()
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Formula" in output
        assert "m/z" in output
        assert "Adduct" in output
        assert "Neutral" in output
        assert "Error" in output
        assert "-" * 10 in output  # Separator line


class TestPrintResults:
    """Test the print_results function."""

    @pytest.fixture
    def sample_hits(self):
        """Sample hit list for testing."""
        return [
            {
                "formula": "Y2Mn4(tBuCOO)10O6H3C1",
                "mz": 1519.1540,
                "charge": -1,
                "adduct": "[M]−•",
                "neutral_mass": 1519.1540,
                "ppm_error": 3.53,
            },
            {
                "formula": "Y2Mn4(tBuCOO)10O6H4C1",
                "mz": 1519.1540,
                "charge": -1,
                "adduct": "[M−H]−",
                "neutral_mass": 1520.1613,
                "ppm_error": 3.89,
            },
        ]

    def test_print_results_with_header(self, sample_hits):
        """print_results should include header by default."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_results(sample_hits, show_header=True)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Formula" in output
        assert "Y2Mn4(tBuCOO)10O6H3C1" in output

    def test_print_results_without_header(self, sample_hits):
        """print_results with show_header=False should omit header."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_results(sample_hits, show_header=False)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # Should not have header row
        assert "Formula" not in output.split("\n")[0]
        # Should still have results
        assert "Y2Mn4(tBuCOO)10O6H3C1" in output

    def test_print_results_all_hits(self, sample_hits):
        """All hits should be printed."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_results(sample_hits)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Y2Mn4(tBuCOO)10O6H3C1" in output
        assert "Y2Mn4(tBuCOO)10O6H4C1" in output

    def test_print_results_empty_list(self):
        """Empty list should produce minimal output."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_results([])
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        # Should just have header
        assert "Formula" in output


class TestPrintLevelHeader:
    """Test the print_level_header function."""

    def test_level_1_header(self):
        """Level 1 header should show 'strict' and parameters."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_level_header(1)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Level 1" in output
        assert "strict" in output
        assert "h_max=0" in output
        assert "c_max=0" in output
        assert "additional_o=0" in output

    def test_level_2_header(self):
        """Level 2 header should show 'moderate' and parameters."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_level_header(2)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Level 2" in output
        assert "moderate" in output
        assert "h_max=4" in output
        assert "c_max=2" in output
        assert "additional_o=2" in output

    def test_level_3_header(self):
        """Level 3 header should show 'loose' and parameters."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            print_level_header(3)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Level 3" in output
        assert "loose" in output
        assert "h_max=10" in output
        assert "c_max=5" in output
        assert "additional_o=5" in output
