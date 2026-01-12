"""Golden file tests - validate search results against known expected outputs.

These tests read expected outputs from JSON fixtures and verify that the
search functions produce matching results. This ensures reproducibility
and catches regressions in the formula enumeration algorithm.
"""

import json
import pytest
from pathlib import Path
from typing import Any

from formula_search import search_mz_negative, search_mz_positive


# Path to fixtures directory
FIXTURES_DIR = Path(__file__).parent / "fixtures"

# Mode configurations: (mode_name, search_function, metal)
MODES = [
    ("y_negative", search_mz_negative, "Y"),
    ("y_positive", search_mz_positive, "Y"),
    ("la_negative", search_mz_negative, "La"),
    ("la_positive", search_mz_positive, "La"),
]


def load_expected_outputs() -> dict:
    """Load expected outputs from JSON fixture file."""
    fixture_path = FIXTURES_DIR / "expected_outputs.json"
    with open(fixture_path) as f:
        return json.load(f)


def get_peak_by_mode(expected_data: dict, mode: str) -> dict | None:
    """Get peak data for a specific mode from fixtures."""
    for peak in expected_data["peaks"]:
        if peak.get("mode") == mode and peak["mz"] == 1519.154:
            return peak
    return None


def normalize_result(result: dict[str, Any]) -> dict[str, Any]:
    """Normalize a result dict for comparison.

    Rounds floating point values and converts counts to list for JSON compatibility.
    """
    return {
        "formula": result["formula"],
        "mz": result["mz"],
        "charge": result["charge"],
        "adduct": result["adduct"],
        "neutral_mass": round(result["neutral_mass"], 6),
        "ppm_error": round(result["ppm_error"], 4),
        "counts": list(result["counts"]),
    }


class TestGoldenOutputs:
    """Test search results against golden (expected) outputs from fixtures."""

    @pytest.fixture
    def expected_data(self) -> dict:
        """Load the expected outputs fixture."""
        return load_expected_outputs()

    def test_fixture_file_exists(self):
        """Verify the fixture file exists and is valid JSON."""
        fixture_path = FIXTURES_DIR / "expected_outputs.json"
        assert fixture_path.exists(), f"Fixture file not found: {fixture_path}"

        with open(fixture_path) as f:
            data = json.load(f)

        assert "peaks" in data
        assert len(data["peaks"]) > 0

    def test_fixture_structure(self, expected_data):
        """Verify fixture has expected structure."""
        for peak in expected_data["peaks"]:
            assert "mz" in peak
            assert "ppm" in peak
            assert "levels" in peak
            assert "mode" in peak
            assert "metal" in peak
            assert "polarity" in peak
            assert all(level in peak["levels"] for level in ["1", "2", "3"])

    def test_all_modes_present(self, expected_data):
        """Verify all 4 modes have fixture data."""
        modes_in_fixture = {p["mode"] for p in expected_data["peaks"]}
        expected_modes = {"y_negative", "y_positive", "la_negative", "la_positive"}
        assert expected_modes.issubset(modes_in_fixture), (
            f"Missing modes: {expected_modes - modes_in_fixture}"
        )

    @pytest.mark.parametrize("level", [1, 2, 3])
    @pytest.mark.parametrize(
        "mode,search_fn,metal",
        MODES,
        ids=["y_negative", "y_positive", "la_negative", "la_positive"],
    )
    def test_peak_1519_all_modes_all_levels(
        self, expected_data, mode, search_fn, metal, level
    ):
        """Test m/z=1519.154 search matches expected output for each mode and level."""
        peak_data = get_peak_by_mode(expected_data, mode)
        assert peak_data is not None, f"No fixture data for mode {mode}"

        mz = peak_data["mz"]
        ppm = peak_data["ppm"]
        expected = peak_data["levels"][str(level)]

        # Run actual search
        actual_results = search_fn(mz, ppm=ppm, coarseness=level, metal=metal)
        actual_normalized = [normalize_result(r) for r in actual_results]

        # Compare counts
        assert len(actual_normalized) == len(expected), (
            f"{mode} level {level}: Expected {len(expected)} results, "
            f"got {len(actual_normalized)}"
        )

        # Compare each result
        for i, (actual, exp) in enumerate(zip(actual_normalized, expected)):
            assert actual["formula"] == exp["formula"], (
                f"{mode} level {level}, result {i}: formula mismatch"
            )
            assert actual["charge"] == exp["charge"], (
                f"{mode} level {level}, result {i}: charge mismatch"
            )
            assert actual["adduct"] == exp["adduct"], (
                f"{mode} level {level}, result {i}: adduct mismatch"
            )
            assert abs(actual["ppm_error"] - exp["ppm_error"]) < 0.001, (
                f"{mode} level {level}, result {i}: ppm_error mismatch "
                f"(got {actual['ppm_error']}, expected {exp['ppm_error']})"
            )
            assert actual["counts"] == exp["counts"], (
                f"{mode} level {level}, result {i}: counts mismatch"
            )


class TestGoldenFormulasYNegative:
    """Test that specific expected formulas are found for Y negative mode."""

    @pytest.fixture
    def expected_data(self) -> dict:
        """Load the expected outputs fixture."""
        return load_expected_outputs()

    def test_y2mn4_radical_anion(self, expected_data):
        """Verify Y2Mn4(tBuCOO)10O6H3C1 [M]-• is found correctly."""
        peak = get_peak_by_mode(expected_data, "y_negative")
        results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=2)

        match = next(
            (
                r
                for r in results
                if r["formula"] == "Y2Mn4(tBuCOO)10O6H3C1" and r["adduct"] == "[M]−•"
            ),
            None,
        )

        assert match is not None, "Y2Mn4(tBuCOO)10O6H3C1 [M]-• not found"
        assert match["charge"] == -1
        assert abs(match["ppm_error"] - 3.531) < 0.01
        assert match["counts"] == (2, 4, 10, 6, 3, 1)

    def test_y2mn4_deprotonated(self, expected_data):
        """Verify Y2Mn4(tBuCOO)10O6H4C1 [M-H]- is found correctly."""
        peak = get_peak_by_mode(expected_data, "y_negative")
        results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=2)

        match = next(
            (
                r
                for r in results
                if r["formula"] == "Y2Mn4(tBuCOO)10O6H4C1" and r["adduct"] == "[M−H]−"
            ),
            None,
        )

        assert match is not None, "Y2Mn4(tBuCOO)10O6H4C1 [M-H]- not found"
        assert match["charge"] == -1
        assert abs(match["ppm_error"] - 3.8895) < 0.01

    def test_y2mn5_chloride_adduct(self, expected_data):
        """Verify Y2Mn5(tBuCOO)10O4H3C2 [M+Cl]- is found correctly."""
        peak = get_peak_by_mode(expected_data, "y_negative")
        results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=2)

        match = next(
            (
                r
                for r in results
                if r["formula"] == "Y2Mn5(tBuCOO)10O4H3C2" and r["adduct"] == "[M+Cl]−"
            ),
            None,
        )

        assert match is not None, "Y2Mn5(tBuCOO)10O4H3C2 [M+Cl]- not found"
        assert abs(match["ppm_error"] - (-9.8279)) < 0.01


class TestGoldenFormulasYPositive:
    """Test that specific expected formulas are found for Y positive mode."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    def test_y1mn5_potassium_adduct(self, expected_data):
        """Verify Y1Mn5(tBuCOO)10O5H2C2 [M+K]+ is found correctly."""
        peak = get_peak_by_mode(expected_data, "y_positive")
        results = search_mz_positive(peak["mz"], ppm=peak["ppm"], coarseness=2, metal="Y")

        match = next(
            (
                r
                for r in results
                if r["formula"] == "Y1Mn5(tBuCOO)10O5H2C2" and r["adduct"] == "[M+K]+"
            ),
            None,
        )

        assert match is not None, "Y1Mn5(tBuCOO)10O5H2C2 [M+K]+ not found"
        assert match["charge"] == 1
        assert abs(match["ppm_error"] - (-0.9898)) < 0.01
        assert match["counts"] == (1, 5, 10, 5, 2, 2)

    def test_y2mn4_protonated(self, expected_data):
        """Verify Y2Mn4(tBuCOO)10O6H2C1 [M+H]+ is found correctly."""
        peak = get_peak_by_mode(expected_data, "y_positive")
        results = search_mz_positive(peak["mz"], ppm=peak["ppm"], coarseness=2, metal="Y")

        match = next(
            (
                r
                for r in results
                if r["formula"] == "Y2Mn4(tBuCOO)10O6H2C1" and r["adduct"] == "[M+H]+"
            ),
            None,
        )

        assert match is not None, "Y2Mn4(tBuCOO)10O6H2C1 [M+H]+ not found"
        assert match["charge"] == 1
        assert abs(match["ppm_error"] - 3.172) < 0.01

    def test_y2mn4_radical_cation(self, expected_data):
        """Verify Y2Mn4(tBuCOO)10O6H3C1 [M]+• is found correctly."""
        peak = get_peak_by_mode(expected_data, "y_positive")
        results = search_mz_positive(peak["mz"], ppm=peak["ppm"], coarseness=2, metal="Y")

        match = next(
            (
                r
                for r in results
                if r["formula"] == "Y2Mn4(tBuCOO)10O6H3C1" and r["adduct"] == "[M]+•"
            ),
            None,
        )

        assert match is not None, "Y2Mn4(tBuCOO)10O6H3C1 [M]+• not found"
        assert match["charge"] == 1
        assert abs(match["ppm_error"] - 3.531) < 0.01


class TestGoldenFormulasLaNegative:
    """Test that specific expected formulas are found for La negative mode."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    def test_la2mn4_deprotonated(self, expected_data):
        """Verify La2Mn4(tBuCOO)10O0H0C1 [M-H]- is found correctly."""
        peak = get_peak_by_mode(expected_data, "la_negative")
        results = search_mz_negative(
            peak["mz"], ppm=peak["ppm"], coarseness=2, metal="La"
        )

        match = next(
            (
                r
                for r in results
                if r["formula"] == "La2Mn4(tBuCOO)10O0H0C1" and r["adduct"] == "[M−H]−"
            ),
            None,
        )

        assert match is not None, "La2Mn4(tBuCOO)10O0H0C1 [M-H]- not found"
        assert match["charge"] == -1
        assert abs(match["ppm_error"] - 2.8844) < 0.01
        assert match["counts"] == (2, 4, 10, 0, 0, 1)

    def test_la1mn5_deprotonated(self, expected_data):
        """Verify La1Mn5(tBuCOO)10O6H0C0 [M-H]- is found correctly."""
        peak = get_peak_by_mode(expected_data, "la_negative")
        results = search_mz_negative(
            peak["mz"], ppm=peak["ppm"], coarseness=2, metal="La"
        )

        match = next(
            (
                r
                for r in results
                if r["formula"] == "La1Mn5(tBuCOO)10O6H0C0" and r["adduct"] == "[M−H]−"
            ),
            None,
        )

        assert match is not None, "La1Mn5(tBuCOO)10O6H0C0 [M-H]- not found"
        assert match["charge"] == -1
        assert abs(match["ppm_error"] - 4.2407) < 0.01


class TestGoldenFormulasLaPositive:
    """Test that specific expected formulas are found for La positive mode."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    def test_la1mn5_potassium_adduct(self, expected_data):
        """Verify La1Mn5(tBuCOO)10O2H0C2 [M+K]+ is found correctly."""
        peak = get_peak_by_mode(expected_data, "la_positive")
        results = search_mz_positive(
            peak["mz"], ppm=peak["ppm"], coarseness=2, metal="La"
        )

        match = next(
            (
                r
                for r in results
                if r["formula"] == "La1Mn5(tBuCOO)10O2H0C2" and r["adduct"] == "[M+K]+"
            ),
            None,
        )

        assert match is not None, "La1Mn5(tBuCOO)10O2H0C2 [M+K]+ not found"
        assert match["charge"] == 1
        assert abs(match["ppm_error"] - (-1.5059)) < 0.01
        assert match["counts"] == (1, 5, 10, 2, 0, 2)

    def test_la2mn4_sodium_adduct_level3(self, expected_data):
        """Verify La2Mn4(tBuCOO)9O2H9C4 [M+Na]+ is found at level 3."""
        peak = get_peak_by_mode(expected_data, "la_positive")
        results = search_mz_positive(
            peak["mz"], ppm=peak["ppm"], coarseness=3, metal="La"
        )

        match = next(
            (
                r
                for r in results
                if r["formula"] == "La2Mn4(tBuCOO)9O2H9C4" and r["adduct"] == "[M+Na]+"
            ),
            None,
        )

        assert match is not None, "La2Mn4(tBuCOO)9O2H9C4 [M+Na]+ not found at level 3"
        assert match["charge"] == 1
        assert abs(match["ppm_error"] - 0.9562) < 0.01


class TestGoldenResultOrder:
    """Test that results are returned in the expected order (by abs ppm error)."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    @pytest.mark.parametrize(
        "mode,search_fn,metal",
        MODES,
        ids=["y_negative", "y_positive", "la_negative", "la_positive"],
    )
    def test_results_sorted_by_ppm(self, expected_data, mode, search_fn, metal):
        """Verify results are sorted by absolute ppm error for all modes."""
        peak = get_peak_by_mode(expected_data, mode)

        for level in [1, 2, 3]:
            results = search_fn(peak["mz"], ppm=peak["ppm"], coarseness=level, metal=metal)

            if len(results) > 1:
                ppm_errors = [abs(r["ppm_error"]) for r in results]
                assert ppm_errors == sorted(ppm_errors), (
                    f"{mode} level {level}: Results not sorted by ppm error"
                )

    @pytest.mark.parametrize(
        "mode,search_fn,metal",
        MODES,
        ids=["y_negative", "y_positive", "la_negative", "la_positive"],
    )
    def test_order_matches_fixture(self, expected_data, mode, search_fn, metal):
        """Verify result order matches fixture exactly for all modes."""
        peak = get_peak_by_mode(expected_data, mode)

        for level in [2, 3]:  # Skip level 1 (often empty)
            expected = peak["levels"][str(level)]
            if not expected:
                continue

            results = search_fn(
                peak["mz"], ppm=peak["ppm"], coarseness=level, metal=metal
            )

            expected_formulas = [r["formula"] for r in expected]
            actual_formulas = [r["formula"] for r in results]

            assert actual_formulas == expected_formulas, (
                f"{mode} level {level}: Result order doesn't match fixture"
            )


class TestGoldenAdductDistribution:
    """Test adduct type distribution matches expected outputs."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    def test_y_negative_level2_adduct_counts(self, expected_data):
        """Verify Y negative level 2 has expected adduct distribution."""
        peak = get_peak_by_mode(expected_data, "y_negative")
        results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=2)

        adduct_counts = {}
        for r in results:
            adduct_counts[r["adduct"]] = adduct_counts.get(r["adduct"], 0) + 1

        # From fixture: 1 radical, 1 deprotonated, 1 chloride
        assert adduct_counts.get("[M]−•", 0) == 1
        assert adduct_counts.get("[M−H]−", 0) == 1
        assert adduct_counts.get("[M+Cl]−", 0) == 1

    def test_y_positive_level2_adduct_counts(self, expected_data):
        """Verify Y positive level 2 has expected adduct distribution."""
        peak = get_peak_by_mode(expected_data, "y_positive")
        results = search_mz_positive(peak["mz"], ppm=peak["ppm"], coarseness=2, metal="Y")

        adduct_counts = {}
        for r in results:
            adduct_counts[r["adduct"]] = adduct_counts.get(r["adduct"], 0) + 1

        # From fixture: 1 potassium, 1 protonated, 1 radical cation
        assert adduct_counts.get("[M+K]+", 0) == 1
        assert adduct_counts.get("[M+H]+", 0) == 1
        assert adduct_counts.get("[M]+•", 0) == 1

    def test_la_negative_level2_adduct_counts(self, expected_data):
        """Verify La negative level 2 has expected adduct distribution."""
        peak = get_peak_by_mode(expected_data, "la_negative")
        results = search_mz_negative(
            peak["mz"], ppm=peak["ppm"], coarseness=2, metal="La"
        )

        adduct_counts = {}
        for r in results:
            adduct_counts[r["adduct"]] = adduct_counts.get(r["adduct"], 0) + 1

        # From fixture: 2 deprotonated
        assert adduct_counts.get("[M−H]−", 0) == 2

    def test_la_positive_level2_adduct_counts(self, expected_data):
        """Verify La positive level 2 has expected adduct distribution."""
        peak = get_peak_by_mode(expected_data, "la_positive")
        results = search_mz_positive(
            peak["mz"], ppm=peak["ppm"], coarseness=2, metal="La"
        )

        adduct_counts = {}
        for r in results:
            adduct_counts[r["adduct"]] = adduct_counts.get(r["adduct"], 0) + 1

        # From fixture: 1 potassium
        assert adduct_counts.get("[M+K]+", 0) == 1


class TestGoldenLevelProgression:
    """Test that higher coarseness levels find more formulas."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    @pytest.mark.parametrize(
        "mode,search_fn,metal",
        MODES,
        ids=["y_negative", "y_positive", "la_negative", "la_positive"],
    )
    def test_level3_has_at_least_level2_results(
        self, expected_data, mode, search_fn, metal
    ):
        """Verify level 3 finds at least as many formulas as level 2."""
        peak = get_peak_by_mode(expected_data, mode)

        results_2 = search_fn(
            peak["mz"], ppm=peak["ppm"], coarseness=2, metal=metal
        )
        results_3 = search_fn(
            peak["mz"], ppm=peak["ppm"], coarseness=3, metal=metal
        )

        assert len(results_3) >= len(results_2), (
            f"{mode}: Level 3 should find at least as many results as level 2"
        )

    def test_y_negative_level3_unique_formulas(self, expected_data):
        """Verify Y negative level 3 finds formulas not present in level 2."""
        peak = get_peak_by_mode(expected_data, "y_negative")

        results_2 = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=2)
        results_3 = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=3)

        formulas_2 = {r["formula"] for r in results_2}
        formulas_3 = {r["formula"] for r in results_3}

        # Level 3 specific formulas from fixture
        level3_unique = {
            "Y2Mn5(tBuCOO)10O3H9C0",
            "Y2Mn5(tBuCOO)10O3H8C0",
            "Y1Mn5(tBuCOO)10O9H1C0",
            "Y1Mn5(tBuCOO)10O9H2C0",
            "Y2Mn5(tBuCOO)10O0H8C4",
            "Y2Mn5(tBuCOO)10O0H9C4",
        }

        for formula in level3_unique:
            assert formula not in formulas_2, f"{formula} should not be in level 2"
            assert formula in formulas_3, f"{formula} should be in level 3"


def generate_golden_file(output_path: Path, peaks: list[dict]) -> None:
    """Utility to regenerate the golden file with new expected outputs.

    Usage:
        from test_golden import generate_golden_file
        generate_golden_file(
            Path("tests/fixtures/expected_outputs.json"),
            [
                {"mz": 1519.154, "ppm": 10, "mode": "y_negative", "metal": "Y",
                 "polarity": "negative", "description": "Example peak"},
                {"mz": 1519.154, "ppm": 10, "mode": "y_positive", "metal": "Y",
                 "polarity": "positive", "description": "Example peak"},
            ]
        )
    """
    fixture = {"peaks": []}

    for peak_config in peaks:
        mz = peak_config["mz"]
        ppm = peak_config["ppm"]
        mode = peak_config["mode"]
        metal = peak_config["metal"]
        polarity = peak_config["polarity"]
        description = peak_config.get("description", f"Peak at m/z={mz} ({mode})")

        # Select search function based on polarity
        search_fn = search_mz_negative if polarity == "negative" else search_mz_positive

        peak_data = {
            "description": description,
            "mz": mz,
            "ppm": ppm,
            "mode": mode,
            "metal": metal,
            "polarity": polarity,
            "levels": {},
        }

        for level in [1, 2, 3]:
            results = search_fn(mz, ppm=ppm, coarseness=level, metal=metal)
            peak_data["levels"][str(level)] = [normalize_result(r) for r in results]

        fixture["peaks"].append(peak_data)

    with open(output_path, "w") as f:
        json.dump(fixture, f, indent=2)

    print(f"Generated golden file: {output_path}")
