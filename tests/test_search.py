"""Tests for formula_search.search module."""

import pytest
from formula_search.search import search_mz_negative, search_mz_positive
from formula_search.constants import ADDUCTS_NEG, ADDUCTS_POS, PROTON, CL35, NA23, K39, NH4


class TestSearchMzNegative:
    """Test the search_mz_negative function."""

    def test_returns_list(self):
        """Function should return a list."""
        results = search_mz_negative(500.0)
        assert isinstance(results, list)

    def test_result_structure(self):
        """Each result should have all required keys."""
        results = search_mz_negative(500.0, ppm=50, coarseness=3)
        if results:
            result = results[0]
            required_keys = [
                "mz",
                "charge",
                "adduct",
                "neutral_mass",
                "formula",
                "mass",
                "ppm_error",
                "counts",
            ]
            for key in required_keys:
                assert key in result, f"Missing key: {key}"

    def test_mz_preserved(self):
        """Input m/z should be preserved in results."""
        mz = 1519.154
        results = search_mz_negative(mz, ppm=10)
        for result in results:
            assert result["mz"] == mz

    def test_charge_in_results(self):
        """Results should have the requested charge state."""
        results = search_mz_negative(500.0, ppm=50, charges=(-1,), coarseness=3)
        for result in results:
            assert result["charge"] == -1

    def test_multiple_charges(self):
        """Should search multiple charge states when requested."""
        results = search_mz_negative(500.0, ppm=100, charges=(-1, -2), coarseness=3)
        charges = set(r["charge"] for r in results)
        # Might not find both, but should include at least requested charges
        assert charges.issubset({-1, -2})

    def test_adduct_in_results(self):
        """Results should include adduct information."""
        results = search_mz_negative(500.0, ppm=50, coarseness=3)
        for result in results:
            assert result["adduct"] in ADDUCTS_NEG

    def test_custom_adducts(self):
        """Should use custom adducts when provided."""
        custom = {"[Test]−": 10.0}
        results = search_mz_negative(500.0, ppm=50, adducts=custom, coarseness=3)
        for result in results:
            assert result["adduct"] == "[Test]−"

    def test_neutral_mass_calculation_deprotonated(self):
        """Neutral mass should be correctly calculated for [M-H]-."""
        mz = 500.0
        results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M−H]−": PROTON}, coarseness=3
        )
        for result in results:
            if result["adduct"] == "[M−H]−":
                expected_neutral = mz + PROTON
                assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_neutral_mass_calculation_radical(self):
        """Neutral mass should be correctly calculated for [M]-•."""
        mz = 500.0
        results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M]−•": 0.0}, coarseness=3
        )
        for result in results:
            expected_neutral = mz  # No mass adjustment
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_neutral_mass_calculation_chloride(self):
        """Neutral mass should be correctly calculated for [M+Cl]-."""
        mz = 500.0
        results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M+Cl]−": CL35}, coarseness=3
        )
        for result in results:
            expected_neutral = mz + CL35
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_results_sorted_by_ppm(self):
        """Results should be sorted by absolute ppm error."""
        results = search_mz_negative(500.0, ppm=50, coarseness=3)
        if len(results) > 1:
            ppm_errors = [abs(r["ppm_error"]) for r in results]
            assert ppm_errors == sorted(ppm_errors)

    def test_max_hits_limit(self):
        """Results should respect max_hits parameter."""
        # Use large ppm to get many hits
        results = search_mz_negative(500.0, ppm=100, max_hits=5, coarseness=3)
        # max_hits is per adduct/charge, so total could be higher
        # but should be limited per combination
        assert len(results) <= 5 * len(ADDUCTS_NEG)

    def test_coarseness_parameter(self):
        """Coarseness should affect results."""
        results_strict = search_mz_negative(500.0, ppm=50, coarseness=1)
        results_loose = search_mz_negative(500.0, ppm=50, coarseness=3)
        # Loose should generally find more
        assert len(results_loose) >= len(results_strict)

    def test_explicit_h_max(self):
        """Explicit h_max should override coarseness default."""
        results = search_mz_negative(500.0, ppm=50, coarseness=1, h_max=10)
        # Should run without error
        assert isinstance(results, list)

    def test_explicit_c_max(self):
        """Explicit c_max should override coarseness default."""
        results = search_mz_negative(500.0, ppm=50, coarseness=1, c_max=5)
        assert isinstance(results, list)

    def test_no_hits_for_very_small_mz(self):
        """Very small m/z shouldn't match valid formulas."""
        results = search_mz_negative(10.0, ppm=1, coarseness=3)
        assert len(results) == 0

    def test_parameter_limits(self):
        """Search parameters should limit results."""
        results = search_mz_negative(
            500.0,
            ppm=50,
            y_max=1,
            mn_max=1,
            tbu_max=3,
            o_max=2,
            coarseness=3,
        )
        for r in results:
            y, mn, k, o, _, _ = r["counts"]
            assert y <= 1
            assert mn <= 1
            assert k <= 3


class TestSearchMzNegativeKnownPeak:
    """Test search with the known peak from user's example: m/z=1519.1540."""

    @pytest.fixture
    def peak_mz(self):
        """The example peak m/z value."""
        return 1519.1540

    def test_finds_results(self, peak_mz):
        """Should find formula matches for the example peak."""
        results = search_mz_negative(peak_mz, ppm=10, coarseness=2)
        assert len(results) > 0

    def test_moderate_level_finds_expected_formulas(self, peak_mz):
        """Level 2 should find specific formulas from the example."""
        results = search_mz_negative(peak_mz, ppm=10, coarseness=2)
        formulas = [r["formula"] for r in results]

        # From the example output (Level 2):
        # Y2Mn4(tBuCOO)10O6H3C1
        # Y2Mn4(tBuCOO)10O6H4C1
        # These should be findable at moderate coarseness
        found_y2mn4 = any("Y2Mn4" in f and "(tBuCOO)10" in f for f in formulas)
        assert found_y2mn4, "Should find Y2Mn4(tBuCOO)10 type formulas"

    def test_loose_level_finds_more_formulas(self, peak_mz):
        """Level 3 should find additional formulas."""
        results_2 = search_mz_negative(peak_mz, ppm=10, coarseness=2)
        results_3 = search_mz_negative(peak_mz, ppm=10, coarseness=3)

        formulas_2 = set(r["formula"] for r in results_2)
        formulas_3 = set(r["formula"] for r in results_3)

        # Level 3 should find at least as many unique formulas
        assert len(formulas_3) >= len(formulas_2)

    def test_finds_different_adducts(self, peak_mz):
        """Should find matches for different adduct types."""
        results = search_mz_negative(peak_mz, ppm=10, coarseness=2)
        adducts_found = set(r["adduct"] for r in results)

        # Should find at least some adduct types
        assert len(adducts_found) >= 1

    def test_ppm_errors_within_tolerance(self, peak_mz):
        """All ppm errors should be within the specified tolerance."""
        ppm = 10
        results = search_mz_negative(peak_mz, ppm=ppm, coarseness=3)
        for r in results:
            assert abs(r["ppm_error"]) <= ppm, (
                f"PPM error {r['ppm_error']} exceeds {ppm}"
            )

    def test_charge_minus_one(self, peak_mz):
        """Default search should use charge -1."""
        results = search_mz_negative(peak_mz, ppm=10)
        for r in results:
            assert r["charge"] == -1


class TestSearchMzPositive:
    """Test the search_mz_positive function."""

    def test_returns_list(self):
        """Function should return a list."""
        results = search_mz_positive(500.0)
        assert isinstance(results, list)

    def test_result_structure(self):
        """Each result should have all required keys."""
        results = search_mz_positive(500.0, ppm=50, coarseness=3)
        if results:
            result = results[0]
            required_keys = [
                "mz",
                "charge",
                "adduct",
                "neutral_mass",
                "formula",
                "mass",
                "ppm_error",
                "counts",
            ]
            for key in required_keys:
                assert key in result, f"Missing key: {key}"

    def test_mz_preserved(self):
        """Input m/z should be preserved in results."""
        mz = 500.0
        results = search_mz_positive(mz, ppm=50, coarseness=3)
        for result in results:
            assert result["mz"] == mz

    def test_charge_in_results(self):
        """Results should have the requested charge state."""
        results = search_mz_positive(500.0, ppm=50, charges=(1,), coarseness=3)
        for result in results:
            assert result["charge"] == 1

    def test_multiple_charges(self):
        """Should search multiple charge states when requested."""
        results = search_mz_positive(500.0, ppm=100, charges=(1, 2), coarseness=3)
        charges = set(r["charge"] for r in results)
        # Might not find both, but should include at least requested charges
        assert charges.issubset({1, 2})

    def test_adduct_in_results(self):
        """Results should include adduct information."""
        results = search_mz_positive(500.0, ppm=50, coarseness=3)
        for result in results:
            assert result["adduct"] in ADDUCTS_POS

    def test_custom_adducts(self):
        """Should use custom adducts when provided."""
        custom = {"[Test]+": 10.0}
        results = search_mz_positive(500.0, ppm=50, adducts=custom, coarseness=3)
        for result in results:
            assert result["adduct"] == "[Test]+"

    def test_neutral_mass_calculation_protonated(self):
        """Neutral mass should be correctly calculated for [M+H]+."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+H]+": PROTON}, coarseness=3
        )
        for result in results:
            if result["adduct"] == "[M+H]+":
                expected_neutral = mz - PROTON
                assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_neutral_mass_calculation_radical(self):
        """Neutral mass should be correctly calculated for [M]+•."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M]+•": 0.0}, coarseness=3
        )
        for result in results:
            expected_neutral = mz  # No mass adjustment
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_neutral_mass_calculation_sodiated(self):
        """Neutral mass should be correctly calculated for [M+Na]+."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+Na]+": NA23}, coarseness=3
        )
        for result in results:
            expected_neutral = mz - NA23
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_neutral_mass_calculation_potassiated(self):
        """Neutral mass should be correctly calculated for [M+K]+."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+K]+": K39}, coarseness=3
        )
        for result in results:
            expected_neutral = mz - K39
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_neutral_mass_calculation_ammoniated(self):
        """Neutral mass should be correctly calculated for [M+NH4]+."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+NH4]+": NH4}, coarseness=3
        )
        for result in results:
            expected_neutral = mz - NH4
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_results_sorted_by_ppm(self):
        """Results should be sorted by absolute ppm error."""
        results = search_mz_positive(500.0, ppm=50, coarseness=3)
        if len(results) > 1:
            ppm_errors = [abs(r["ppm_error"]) for r in results]
            assert ppm_errors == sorted(ppm_errors)

    def test_max_hits_limit(self):
        """Results should respect max_hits parameter."""
        # Use large ppm to get many hits
        results = search_mz_positive(500.0, ppm=100, max_hits=5, coarseness=3)
        # max_hits is per adduct/charge, so total could be higher
        # but should be limited per combination
        assert len(results) <= 5 * len(ADDUCTS_POS)

    def test_coarseness_parameter(self):
        """Coarseness should affect results."""
        results_strict = search_mz_positive(500.0, ppm=50, coarseness=1)
        results_loose = search_mz_positive(500.0, ppm=50, coarseness=3)
        # Loose should generally find more
        assert len(results_loose) >= len(results_strict)

    def test_explicit_h_max(self):
        """Explicit h_max should override coarseness default."""
        results = search_mz_positive(500.0, ppm=50, coarseness=1, h_max=10)
        # Should run without error
        assert isinstance(results, list)

    def test_explicit_c_max(self):
        """Explicit c_max should override coarseness default."""
        results = search_mz_positive(500.0, ppm=50, coarseness=1, c_max=5)
        assert isinstance(results, list)

    def test_no_hits_for_very_small_mz(self):
        """Very small m/z shouldn't match valid formulas."""
        results = search_mz_positive(10.0, ppm=1, coarseness=3)
        assert len(results) == 0

    def test_parameter_limits(self):
        """Search parameters should limit results."""
        results = search_mz_positive(
            500.0,
            ppm=50,
            y_max=1,
            mn_max=1,
            tbu_max=3,
            o_max=2,
            coarseness=3,
        )
        for r in results:
            y, mn, k, o, _, _ = r["counts"]
            assert y <= 1
            assert mn <= 1
            assert k <= 3

    def test_charge_plus_one(self):
        """Default search should use charge +1."""
        results = search_mz_positive(500.0, ppm=50, coarseness=3)
        for r in results:
            assert r["charge"] == 1

    def test_ppm_errors_within_tolerance(self):
        """All ppm errors should be within the specified tolerance."""
        ppm = 10
        results = search_mz_positive(500.0, ppm=ppm, coarseness=3)
        for r in results:
            assert abs(r["ppm_error"]) <= ppm, (
                f"PPM error {r['ppm_error']} exceeds {ppm}"
            )


class TestSearchLanthanumSupport:
    """Test La (lanthanum) based complex searches."""

    def test_negative_mode_la_metal(self):
        """Negative mode search should support La metal."""
        results = search_mz_negative(500.0, ppm=50, coarseness=3, metal="La")
        assert isinstance(results, list)
        for result in results:
            assert "La" in result["formula"]
            assert result.get("metal") == "La"

    def test_positive_mode_la_metal(self):
        """Positive mode search should support La metal."""
        results = search_mz_positive(500.0, ppm=50, coarseness=3, metal="La")
        assert isinstance(results, list)
        for result in results:
            assert "La" in result["formula"]
            assert result.get("metal") == "La"

    def test_default_metal_is_y(self):
        """Default metal should be Y."""
        results = search_mz_negative(500.0, ppm=50, coarseness=3)
        for result in results:
            assert "Y" in result["formula"]

    def test_la_vs_y_different_results(self):
        """La and Y searches for same m/z should give different results."""
        y_results = search_mz_negative(500.0, ppm=50, coarseness=3, metal="Y")
        la_results = search_mz_negative(500.0, ppm=50, coarseness=3, metal="La")

        y_formulas = set(r["formula"] for r in y_results)
        la_formulas = set(r["formula"] for r in la_results)

        # Formulas should be different (different metal prefixes)
        assert y_formulas.isdisjoint(la_formulas) or len(y_formulas) == 0


class TestSearchMzPositiveVsNegative:
    """Test differences between positive and negative mode searches."""

    def test_different_neutral_mass_for_same_mz(self):
        """Positive and negative mode should calculate different neutral masses."""
        mz = 500.0
        pos_results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+H]+": PROTON}, coarseness=3
        )
        neg_results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M−H]−": PROTON}, coarseness=3
        )

        if pos_results and neg_results:
            # [M+H]+: neutral = mz - proton
            # [M-H]-: neutral = mz + proton
            pos_neutral = pos_results[0]["neutral_mass"]
            neg_neutral = neg_results[0]["neutral_mass"]
            # Difference should be approximately 2 * proton mass
            assert abs((neg_neutral - pos_neutral) - 2 * PROTON) < 1e-6

    def test_positive_mode_default_adducts(self):
        """Positive mode should use ADDUCTS_POS by default."""
        results = search_mz_positive(500.0, ppm=50, coarseness=3)
        adducts_used = set(r["adduct"] for r in results)
        assert adducts_used.issubset(set(ADDUCTS_POS.keys()))

    def test_negative_mode_default_adducts(self):
        """Negative mode should use ADDUCTS_NEG by default."""
        results = search_mz_negative(500.0, ppm=50, coarseness=3)
        adducts_used = set(r["adduct"] for r in results)
        assert adducts_used.issubset(set(ADDUCTS_NEG.keys()))
