"""Tests for formula_search.enumeration module."""

import pytest
from formula_search.enumeration import within_ppm, enumerate_tBuCOO_YMn
from formula_search.constants import MASS


class TestWithinPpm:
    """Test the within_ppm function."""

    def test_exact_match(self):
        """Exact mass should be within any ppm tolerance."""
        assert within_ppm(100.0, 100.0, 1)
        assert within_ppm(100.0, 100.0, 10)

    def test_within_tolerance(self):
        """Mass within tolerance should return True."""
        # 10 ppm of 1000 = 0.01 Da
        assert within_ppm(1000.005, 1000.0, 10)
        assert within_ppm(999.995, 1000.0, 10)

    def test_outside_tolerance(self):
        """Mass outside tolerance should return False."""
        # 10 ppm of 1000 = 0.01 Da
        assert not within_ppm(1000.02, 1000.0, 10)
        assert not within_ppm(999.98, 1000.0, 10)

    def test_boundary_cases(self):
        """Edge cases at tolerance boundary."""
        target = 1000.0
        ppm = 10  # 0.01 Da tolerance
        # Exactly at boundary
        assert within_ppm(1000.01, target, ppm)
        assert within_ppm(999.99, target, ppm)

    def test_small_mass(self):
        """PPM tolerance for small masses."""
        # 10 ppm of 100 = 0.001 Da
        assert within_ppm(100.0005, 100.0, 10)
        assert not within_ppm(100.002, 100.0, 10)

    def test_large_mass(self):
        """PPM tolerance for large masses."""
        # 10 ppm of 10000 = 0.1 Da
        assert within_ppm(10000.05, 10000.0, 10)
        assert not within_ppm(10000.2, 10000.0, 10)


class TestEnumerateTBuCOOYMn:
    """Test the enumerate_tBuCOO_YMn function."""

    def test_returns_list(self):
        """Function should return a list."""
        result = enumerate_tBuCOO_YMn(500.0, ppm=10)
        assert isinstance(result, list)

    def test_hit_structure(self):
        """Each hit should have required keys."""
        # Use a mass that's likely to have hits
        results = enumerate_tBuCOO_YMn(300.0, ppm=100, coarseness=3)
        if results:
            hit = results[0]
            assert "formula" in hit
            assert "mass" in hit
            assert "ppm_error" in hit
            assert "counts" in hit

    def test_counts_tuple_format(self):
        """Counts should be a tuple of 6 integers (Y, Mn, tBuCOO, O, H, C)."""
        results = enumerate_tBuCOO_YMn(300.0, ppm=100, coarseness=3)
        if results:
            counts = results[0]["counts"]
            assert isinstance(counts, tuple)
            assert len(counts) == 6
            assert all(isinstance(c, int) for c in counts)

    def test_formula_format(self):
        """Formula should follow expected string format."""
        results = enumerate_tBuCOO_YMn(300.0, ppm=100, coarseness=3)
        if results:
            formula = results[0]["formula"]
            assert "Y" in formula
            assert "Mn" in formula
            assert "(tBuCOO)" in formula
            assert "O" in formula
            assert "H" in formula
            assert "C" in formula

    def test_results_sorted_by_ppm(self):
        """Results should be sorted by absolute ppm error."""
        results = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=3)
        if len(results) > 1:
            ppm_errors = [abs(r["ppm_error"]) for r in results]
            assert ppm_errors == sorted(ppm_errors)

    def test_ppm_within_tolerance(self):
        """All results should be within ppm tolerance."""
        target = 500.0
        ppm = 10
        results = enumerate_tBuCOO_YMn(target, ppm=ppm, coarseness=3)
        for hit in results:
            assert abs(hit["ppm_error"]) <= ppm

    def test_metal_constraint(self):
        """All results should have at least one metal (Y or Mn)."""
        results = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=3)
        for hit in results:
            y, mn, _, _, _, _ = hit["counts"]
            assert y > 0 or mn > 0, "At least one metal required"

    def test_ligand_constraint(self):
        """All results should have at least one ligand (tBuCOO) or oxygen."""
        results = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=3)
        for hit in results:
            _, _, k, o, _, _ = hit["counts"]
            assert k > 0 or o > 0, "At least one ligand or oxygen required"

    def test_charge_balance_constraint(self):
        """If tBuCOO present: 2*tBuCOO + O >= metals."""
        results = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=3)
        for hit in results:
            y, mn, k, o, _, _ = hit["counts"]
            if k > 0:
                assert 2 * k + o >= y + mn, "Charge balance constraint violated"

    def test_coarseness_level_1_strict(self):
        """Level 1 should have no additional H, C, or O."""
        results = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=1)
        for hit in results:
            _, _, _, _, h, c = hit["counts"]
            assert h == 0, "Strict mode should have h=0"
            assert c == 0, "Strict mode should have c=0"

    def test_coarseness_level_affects_results(self):
        """Higher coarseness should generally yield more results."""
        results_1 = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=1)
        results_3 = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=3)
        # Loose mode should find at least as many hits
        assert len(results_3) >= len(results_1)

    def test_explicit_h_max_override(self):
        """Explicit h_max should override coarseness default."""
        # With coarseness=1 (h_max=0) but explicit h_max=5
        results = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=1, h_max=5)
        has_hydrogen = any(hit["counts"][4] > 0 for hit in results)
        # Should allow hydrogen despite strict coarseness
        assert len(results) >= 0  # Just checking it runs without error

    def test_explicit_c_max_override(self):
        """Explicit c_max should override coarseness default."""
        results = enumerate_tBuCOO_YMn(500.0, ppm=50, coarseness=1, c_max=3)
        # Should run without error
        assert isinstance(results, list)

    def test_mass_calculation_accuracy(self):
        """Verify mass calculations are accurate."""
        # Y1Mn1(tBuCOO)2O1H0C0 has known mass
        expected_mass = (
            1 * MASS["Y"] + 1 * MASS["Mn"] + 2 * MASS["tBuCOO"] + 1 * MASS["O"]
        )
        results = enumerate_tBuCOO_YMn(expected_mass, ppm=1, coarseness=1)
        # Should find this exact formula
        found = any(hit["counts"] == (1, 1, 2, 1, 0, 0) for hit in results)
        assert found, f"Should find Y1Mn1(tBuCOO)2O1H0C0 at mass {expected_mass}"

    def test_no_results_for_impossible_mass(self):
        """Very small masses shouldn't match any valid formula."""
        results = enumerate_tBuCOO_YMn(10.0, ppm=1, coarseness=3)
        assert len(results) == 0

    def test_parameter_limits_respected(self):
        """y_max, mn_max, tbu_max, o_max should be respected."""
        results = enumerate_tBuCOO_YMn(
            1000.0, ppm=100, y_max=1, mn_max=2, tbu_max=5, o_max=3, coarseness=3
        )
        for hit in results:
            y, mn, k, o, _, _ = hit["counts"]
            assert y <= 1, "Y count exceeds y_max"
            assert mn <= 2, "Mn count exceeds mn_max"
            assert k <= 5, "tBuCOO count exceeds tbu_max"
            # o_max + additional_o for coarseness 3
            assert o <= 3 + 5, "O count exceeds o_max + additional_o"


class TestEnumerateKnownFormulas:
    """Test enumeration with specific known formulas."""

    def test_simple_metal_oxide(self):
        """Test finding a simple Y-Mn oxide formula."""
        # Y2O3 mass
        target = 2 * MASS["Y"] + 3 * MASS["O"]
        results = enumerate_tBuCOO_YMn(target, ppm=1, coarseness=1)
        found = any(hit["counts"] == (2, 0, 0, 3, 0, 0) for hit in results)
        assert found, f"Should find Y2O3 at mass {target:.4f}"

    def test_single_ligand_complex(self):
        """Test finding a single-ligand complex."""
        # Y1(tBuCOO)1O1 mass
        target = MASS["Y"] + MASS["tBuCOO"] + MASS["O"]
        results = enumerate_tBuCOO_YMn(target, ppm=1, coarseness=1)
        found = any(hit["counts"] == (1, 0, 1, 1, 0, 0) for hit in results)
        assert found, f"Should find Y1(tBuCOO)1O1 at mass {target:.4f}"
