"""Integration tests for formula_search package.

These tests verify the complete workflow matches expected outputs
from the user's example with m/z=1519.1540.
"""

import pytest
import io
import sys
from formula_search import (
    search_mz_negative,
    format_hit,
    run_scan_all,
)


class TestIntegrationExample:
    """Integration tests based on the user's example output."""

    @pytest.fixture
    def example_mz(self):
        """The example m/z value from user's output."""
        return 1519.1540

    def test_full_search_workflow(self, example_mz):
        """Test complete search and formatting workflow."""
        # Search
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        
        # Verify results exist
        assert len(results) > 0
        
        # Format each result
        for r in results:
            formatted = format_hit(r)
            assert isinstance(formatted, str)
            assert len(formatted) > 40

    def test_example_output_format_matches(self, example_mz):
        """Verify output format matches the expected pattern."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        
        if results:
            formatted = format_hit(results[0])
            
            # Expected format (from user example):
            # Y2Mn4(tBuCOO)10O6H3C1 mz=1519.1540 charge=-1 adduct=[M]−• neutral_mass=1519.1540 ppm_error=3.53
            
            # Verify key components present
            assert "mz" not in formatted.lower() or "1519" in formatted
            assert results[0]["formula"] in formatted
            assert "ppm" in formatted

    def test_scan_all_workflow(self, example_mz):
        """Test the run_scan_all function produces output."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            run_scan_all(example_mz, ppm=10)
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        
        # Should have level headers
        assert "Level 1" in output or "strict" in output
        assert "Level 2" in output or "moderate" in output
        assert "Level 3" in output or "loose" in output
        
        # Should have some formulas
        assert "Y" in output
        assert "Mn" in output

    def test_incremental_level_results(self, example_mz):
        """Each level should include all previous level results."""
        results_1 = search_mz_negative(example_mz, ppm=10, coarseness=1)
        results_2 = search_mz_negative(example_mz, ppm=10, coarseness=2)
        results_3 = search_mz_negative(example_mz, ppm=10, coarseness=3)
        
        formulas_1 = set(r["formula"] for r in results_1)
        formulas_2 = set(r["formula"] for r in results_2)
        formulas_3 = set(r["formula"] for r in results_3)
        
        # Higher levels should include lower level formulas
        assert formulas_1.issubset(formulas_2)
        assert formulas_2.issubset(formulas_3)


class TestExpectedFormulas:
    """Test that specific expected formulas are found."""

    @pytest.fixture
    def example_mz(self):
        return 1519.1540

    def test_finds_y2mn4_formulas(self, example_mz):
        """Should find Y2Mn4 type formulas from the example."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        formulas = [r["formula"] for r in results]
        
        # From example: Y2Mn4(tBuCOO)10O6H3C1, Y2Mn4(tBuCOO)10O6H4C1
        y2mn4_found = any("Y2" in f and "Mn4" in f for f in formulas)
        assert y2mn4_found, f"Expected Y2Mn4 formulas, got: {formulas[:5]}"

    def test_finds_tBuCOO10_formulas(self, example_mz):
        """Should find formulas with 10 tBuCOO ligands."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        formulas = [r["formula"] for r in results]
        
        tbu10_found = any("(tBuCOO)10" in f for f in formulas)
        assert tbu10_found, f"Expected (tBuCOO)10 formulas, got: {formulas[:5]}"

    def test_radical_anion_adduct_found(self, example_mz):
        """Should find [M]-• radical anion matches."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        adducts = [r["adduct"] for r in results]
        
        assert "[M]−•" in adducts, f"Expected [M]−• adduct, got: {set(adducts)}"

    def test_deprotonated_adduct_found(self, example_mz):
        """Should find [M-H]- deprotonated matches."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        adducts = [r["adduct"] for r in results]
        
        assert "[M−H]−" in adducts, f"Expected [M−H]− adduct, got: {set(adducts)}"

    def test_level3_finds_y1mn5_formulas(self, example_mz):
        """Level 3 should find Y1Mn5 formulas from the example."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=3)
        formulas = [r["formula"] for r in results]
        
        # From example: Y1Mn5(tBuCOO)10O9H1C0
        y1mn5_found = any("Y1" in f and "Mn5" in f for f in formulas)
        assert y1mn5_found, f"Expected Y1Mn5 formulas at level 3"

    def test_level3_finds_y2mn5_formulas(self, example_mz):
        """Level 3 should find Y2Mn5 formulas from the example."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=3)
        formulas = [r["formula"] for r in results]
        
        # From example: Y2Mn5(tBuCOO)10O3H9C0
        y2mn5_found = any("Y2" in f and "Mn5" in f for f in formulas)
        assert y2mn5_found, f"Expected Y2Mn5 formulas at level 3"


class TestPpmErrorValues:
    """Test that ppm errors match expected ranges."""

    @pytest.fixture
    def example_mz(self):
        return 1519.1540

    def test_best_match_ppm_error(self, example_mz):
        """Best match should have low ppm error."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        
        if results:
            best_ppm = abs(results[0]["ppm_error"])
            # From example, best matches are ~3-4 ppm
            assert best_ppm < 10, f"Best ppm error {best_ppm} should be < 10"

    def test_all_results_within_tolerance(self, example_mz):
        """All results should be within specified ppm tolerance."""
        ppm = 10
        results = search_mz_negative(example_mz, ppm=ppm, coarseness=3)
        
        for r in results:
            assert abs(r["ppm_error"]) <= ppm, (
                f"Result {r['formula']} has ppm_error {r['ppm_error']} "
                f"exceeding tolerance {ppm}"
            )


class TestMultipleAdductTypes:
    """Test that multiple adduct types are correctly handled."""

    @pytest.fixture
    def example_mz(self):
        return 1519.1540

    def test_chloride_adduct_produces_different_neutral(self, example_mz):
        """[M+Cl]- should give higher neutral mass than [M]-•."""
        results = search_mz_negative(example_mz, ppm=15, coarseness=2)
        
        radical_results = [r for r in results if r["adduct"] == "[M]−•"]
        chloride_results = [r for r in results if r["adduct"] == "[M+Cl]−"]
        
        if radical_results and chloride_results:
            # Chloride adduct neutral mass should be ~35 Da higher
            radical_neutral = radical_results[0]["neutral_mass"]
            chloride_neutral = chloride_results[0]["neutral_mass"]
            
            assert chloride_neutral > radical_neutral + 30, (
                f"Chloride neutral {chloride_neutral} should be "
                f">35 Da higher than radical {radical_neutral}"
            )

    def test_deprotonated_adduct_neutral_mass(self, example_mz):
        """[M-H]- should give neutral mass = m/z + proton mass."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)
        
        deprotonated = [r for r in results if r["adduct"] == "[M−H]−"]
        radical = [r for r in results if r["adduct"] == "[M]−•"]
        
        if deprotonated and radical:
            # Deprotonated neutral should be ~1 Da higher than radical
            deprot_neutral = deprotonated[0]["neutral_mass"]
            radical_neutral = radical[0]["neutral_mass"]
            
            assert 0.9 < (deprot_neutral - radical_neutral) < 1.1, (
                f"Deprotonated - radical neutral mass difference "
                f"should be ~1 Da (proton)"
            )
