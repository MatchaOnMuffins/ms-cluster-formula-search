"""Tests for formula_search.search module."""

import pytest
from formula_search.search import search_mz_negative, search_mz_positive
from formula_search.constants import ADDUCTS_NEG, ADDUCTS_POS, PROTON, CL35, NA23, K39, NH4


# Parametrized test configurations
SEARCH_MODES = [
    pytest.param(search_mz_negative, ADDUCTS_NEG, -1, id="negative"),
    pytest.param(search_mz_positive, ADDUCTS_POS, 1, id="positive"),
]


class TestSearchBasics:
    """Test basic search function behavior for both modes."""

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_returns_list(self, search_fn, adducts, charge):
        """Function should return a list."""
        results = search_fn(500.0)
        assert isinstance(results, list)

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_result_structure(self, search_fn, adducts, charge):
        """Each result should have all required keys."""
        results = search_fn(500.0, ppm=50, coarseness=3)
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

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_mz_preserved(self, search_fn, adducts, charge):
        """Input m/z should be preserved in results."""
        mz = 500.0
        results = search_fn(mz, ppm=50, coarseness=3)
        for result in results:
            assert result["mz"] == mz

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_charge_in_results(self, search_fn, adducts, charge):
        """Results should have the requested charge state."""
        results = search_fn(500.0, ppm=50, charges=(charge,), coarseness=3)
        for result in results:
            assert result["charge"] == charge

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_adduct_in_results(self, search_fn, adducts, charge):
        """Results should include valid adduct information."""
        results = search_fn(500.0, ppm=50, coarseness=3)
        for result in results:
            assert result["adduct"] in adducts

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_custom_adducts(self, search_fn, adducts, charge):
        """Should use custom adducts when provided."""
        sign = "+" if charge > 0 else "−"
        custom = {f"[Test]{sign}": 10.0}
        results = search_fn(500.0, ppm=50, adducts=custom, coarseness=3)
        for result in results:
            assert result["adduct"] == f"[Test]{sign}"

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_results_sorted_by_ppm(self, search_fn, adducts, charge):
        """Results should be sorted by absolute ppm error."""
        results = search_fn(500.0, ppm=50, coarseness=3)
        if len(results) > 1:
            ppm_errors = [abs(r["ppm_error"]) for r in results]
            assert ppm_errors == sorted(ppm_errors)

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_max_hits_limit(self, search_fn, adducts, charge):
        """Results should respect max_hits parameter."""
        results = search_fn(500.0, ppm=100, max_hits=5, coarseness=3)
        assert len(results) <= 5 * len(adducts)

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_coarseness_parameter(self, search_fn, adducts, charge):
        """Higher coarseness should find at least as many results."""
        results_strict = search_fn(500.0, ppm=50, coarseness=1)
        results_loose = search_fn(500.0, ppm=50, coarseness=3)
        assert len(results_loose) >= len(results_strict)

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_explicit_h_max(self, search_fn, adducts, charge):
        """Explicit h_max should override coarseness default."""
        results = search_fn(500.0, ppm=50, coarseness=1, h_max=10)
        assert isinstance(results, list)

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_explicit_c_max(self, search_fn, adducts, charge):
        """Explicit c_max should override coarseness default."""
        results = search_fn(500.0, ppm=50, coarseness=1, c_max=5)
        assert isinstance(results, list)

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_no_hits_for_very_small_mz(self, search_fn, adducts, charge):
        """Very small m/z shouldn't match valid formulas."""
        results = search_fn(10.0, ppm=1, coarseness=3)
        assert len(results) == 0

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_ppm_errors_within_tolerance(self, search_fn, adducts, charge):
        """All ppm errors should be within the specified tolerance."""
        ppm = 10
        results = search_fn(500.0, ppm=ppm, coarseness=3)
        for r in results:
            assert abs(r["ppm_error"]) <= ppm, (
                f"PPM error {r['ppm_error']} exceeds {ppm}"
            )

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_parameter_limits(self, search_fn, adducts, charge):
        """Search parameters should limit results."""
        results = search_fn(
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


class TestNeutralMassCalculations:
    """Test neutral mass calculations for various adducts."""

    def test_negative_deprotonated(self):
        """[M-H]- neutral mass should be mz + proton."""
        mz = 500.0
        results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M−H]−": PROTON}, coarseness=3
        )
        for result in results:
            expected_neutral = mz + PROTON
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_negative_radical(self):
        """[M]-• neutral mass should equal mz (no adjustment)."""
        mz = 500.0
        results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M]−•": 0.0}, coarseness=3
        )
        for result in results:
            assert abs(result["neutral_mass"] - mz) < 1e-6

    def test_negative_chloride(self):
        """[M+Cl]- neutral mass should be mz + Cl mass."""
        mz = 500.0
        results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M+Cl]−": CL35}, coarseness=3
        )
        for result in results:
            expected_neutral = mz + CL35
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_positive_protonated(self):
        """[M+H]+ neutral mass should be mz - proton."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+H]+": PROTON}, coarseness=3
        )
        for result in results:
            expected_neutral = mz - PROTON
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_positive_radical(self):
        """[M]+• neutral mass should equal mz (no adjustment)."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M]+•": 0.0}, coarseness=3
        )
        for result in results:
            assert abs(result["neutral_mass"] - mz) < 1e-6

    def test_positive_sodiated(self):
        """[M+Na]+ neutral mass should be mz - Na mass."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+Na]+": NA23}, coarseness=3
        )
        for result in results:
            expected_neutral = mz - NA23
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_positive_potassiated(self):
        """[M+K]+ neutral mass should be mz - K mass."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+K]+": K39}, coarseness=3
        )
        for result in results:
            expected_neutral = mz - K39
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_positive_ammoniated(self):
        """[M+NH4]+ neutral mass should be mz - NH4 mass."""
        mz = 500.0
        results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+NH4]+": NH4}, coarseness=3
        )
        for result in results:
            expected_neutral = mz - NH4
            assert abs(result["neutral_mass"] - expected_neutral) < 1e-6

    def test_positive_vs_negative_neutral_mass_difference(self):
        """Same mz with [M+H]+ vs [M-H]- should differ by ~2*proton."""
        mz = 500.0
        pos_results = search_mz_positive(
            mz, ppm=50, charges=(1,), adducts={"[M+H]+": PROTON}, coarseness=3
        )
        neg_results = search_mz_negative(
            mz, ppm=50, charges=(-1,), adducts={"[M−H]−": PROTON}, coarseness=3
        )

        if pos_results and neg_results:
            pos_neutral = pos_results[0]["neutral_mass"]
            neg_neutral = neg_results[0]["neutral_mass"]
            assert abs((neg_neutral - pos_neutral) - 2 * PROTON) < 1e-6


class TestMultipleCharges:
    """Test multiple charge state searches."""

    def test_negative_multiple_charges(self):
        """Should search multiple charge states when requested."""
        results = search_mz_negative(500.0, ppm=100, charges=(-1, -2), coarseness=3)
        charges = set(r["charge"] for r in results)
        assert charges.issubset({-1, -2})

    def test_positive_multiple_charges(self):
        """Should search multiple charge states when requested."""
        results = search_mz_positive(500.0, ppm=100, charges=(1, 2), coarseness=3)
        charges = set(r["charge"] for r in results)
        assert charges.issubset({1, 2})


class TestMetalSupport:
    """Test Y and La metal searches."""

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_default_metal_is_y(self, search_fn, adducts, charge):
        """Default metal should be Y."""
        results = search_fn(500.0, ppm=50, coarseness=3)
        for result in results:
            assert "Y" in result["formula"]

    @pytest.mark.parametrize("search_fn,adducts,charge", SEARCH_MODES)
    def test_la_metal_support(self, search_fn, adducts, charge):
        """Should support La metal searches."""
        results = search_fn(500.0, ppm=50, coarseness=3, metal="La")
        assert isinstance(results, list)
        for result in results:
            assert "La" in result["formula"]
            assert result.get("metal") == "La"

    def test_la_vs_y_different_results(self):
        """La and Y searches for same m/z should give different formulas."""
        y_results = search_mz_negative(500.0, ppm=50, coarseness=3, metal="Y")
        la_results = search_mz_negative(500.0, ppm=50, coarseness=3, metal="La")

        y_formulas = set(r["formula"] for r in y_results)
        la_formulas = set(r["formula"] for r in la_results)

        # Formulas should be different (different metal prefixes)
        assert y_formulas.isdisjoint(la_formulas) or len(y_formulas) == 0
