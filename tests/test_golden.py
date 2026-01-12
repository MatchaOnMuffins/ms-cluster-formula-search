"""Golden file tests - validate search results against known expected outputs.

These tests read expected outputs from JSON fixtures and verify that the
search functions produce matching results. This ensures reproducibility
and catches regressions in the formula enumeration algorithm.
"""

import json
import pytest
from pathlib import Path
from typing import Any

from formula_search import search_mz_negative


# Path to fixtures directory
FIXTURES_DIR = Path(__file__).parent / "fixtures"


def load_expected_outputs() -> dict:
    """Load expected outputs from JSON fixture file."""
    fixture_path = FIXTURES_DIR / "expected_outputs.json"
    with open(fixture_path) as f:
        return json.load(f)


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
            assert all(level in peak["levels"] for level in ["1", "2", "3"])

    @pytest.mark.parametrize("level", [1, 2, 3])
    def test_peak_1519_all_levels(self, expected_data, level):
        """Test m/z=1519.154 search matches expected output at each level."""
        # Find the 1519.154 peak in fixtures
        peak_data = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

        mz = peak_data["mz"]
        ppm = peak_data["ppm"]
        expected = peak_data["levels"][str(level)]

        # Run actual search
        actual_results = search_mz_negative(mz, ppm=ppm, coarseness=level)
        actual_normalized = [normalize_result(r) for r in actual_results]

        # Compare counts
        assert len(actual_normalized) == len(expected), (
            f"Level {level}: Expected {len(expected)} results, got {len(actual_normalized)}"
        )

        # Compare each result
        for i, (actual, exp) in enumerate(zip(actual_normalized, expected)):
            assert actual["formula"] == exp["formula"], (
                f"Level {level}, result {i}: formula mismatch"
            )
            assert actual["charge"] == exp["charge"], (
                f"Level {level}, result {i}: charge mismatch"
            )
            assert actual["adduct"] == exp["adduct"], (
                f"Level {level}, result {i}: adduct mismatch"
            )
            assert abs(actual["ppm_error"] - exp["ppm_error"]) < 0.001, (
                f"Level {level}, result {i}: ppm_error mismatch "
                f"(got {actual['ppm_error']}, expected {exp['ppm_error']})"
            )
            assert actual["counts"] == exp["counts"], (
                f"Level {level}, result {i}: counts mismatch"
            )


class TestGoldenFormulas:
    """Test that specific expected formulas are found with correct properties."""

    @pytest.fixture
    def expected_data(self) -> dict:
        """Load the expected outputs fixture."""
        return load_expected_outputs()

    def test_y2mn4_radical_anion(self, expected_data):
        """Verify Y2Mn4(tBuCOO)10O6H3C1 [M]-• is found correctly."""
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

        results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=2)

        # Find the specific formula
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
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

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
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

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

    def test_level3_unique_formulas(self, expected_data):
        """Verify level 3 finds formulas not present in level 2."""
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

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


class TestGoldenResultOrder:
    """Test that results are returned in the expected order (by abs ppm error)."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    def test_results_sorted_by_ppm(self, expected_data):
        """Verify results are sorted by absolute ppm error."""
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

        for level in [1, 2, 3]:
            results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=level)

            if len(results) > 1:
                ppm_errors = [abs(r["ppm_error"]) for r in results]
                assert ppm_errors == sorted(ppm_errors), (
                    f"Level {level}: Results not sorted by ppm error"
                )

    def test_order_matches_fixture(self, expected_data):
        """Verify result order matches fixture exactly."""
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

        for level in [2, 3]:  # Skip level 1 (empty)
            expected = peak["levels"][str(level)]
            results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=level)

            expected_formulas = [r["formula"] for r in expected]
            actual_formulas = [r["formula"] for r in results]

            assert actual_formulas == expected_formulas, (
                f"Level {level}: Result order doesn't match fixture"
            )


class TestGoldenAdductDistribution:
    """Test adduct type distribution matches expected outputs."""

    @pytest.fixture
    def expected_data(self) -> dict:
        return load_expected_outputs()

    def test_level2_adduct_counts(self, expected_data):
        """Verify level 2 has expected adduct distribution."""
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

        results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=2)

        adduct_counts = {}
        for r in results:
            adduct_counts[r["adduct"]] = adduct_counts.get(r["adduct"], 0) + 1

        # From fixture: 1 radical, 1 deprotonated, 1 chloride
        assert adduct_counts.get("[M]−•", 0) == 1
        assert adduct_counts.get("[M−H]−", 0) == 1
        assert adduct_counts.get("[M+Cl]−", 0) == 1

    def test_level3_adduct_counts(self, expected_data):
        """Verify level 3 has expected adduct distribution."""
        peak = next(
            p for p in expected_data["peaks"] if abs(p["mz"] - 1519.154) < 0.001
        )

        results = search_mz_negative(peak["mz"], ppm=peak["ppm"], coarseness=3)

        adduct_counts = {}
        for r in results:
            adduct_counts[r["adduct"]] = adduct_counts.get(r["adduct"], 0) + 1

        # From fixture: 4 radical, 4 deprotonated, 2 chloride
        assert adduct_counts.get("[M]−•", 0) == 4
        assert adduct_counts.get("[M−H]−", 0) == 4
        assert adduct_counts.get("[M+Cl]−", 0) == 2


def generate_golden_file(output_path: Path, peaks: list[dict]) -> None:
    """Utility to regenerate the golden file with new expected outputs.

    Usage:
        from test_golden import generate_golden_file
        generate_golden_file(
            Path("tests/fixtures/expected_outputs.json"),
            [{"mz": 1519.154, "ppm": 10, "description": "Example peak"}]
        )
    """
    fixture = {"peaks": []}

    for peak_config in peaks:
        mz = peak_config["mz"]
        ppm = peak_config["ppm"]
        description = peak_config.get("description", f"Peak at m/z={mz}")

        peak_data = {"description": description, "mz": mz, "ppm": ppm, "levels": {}}

        for level in [1, 2, 3]:
            results = search_mz_negative(mz, ppm=ppm, coarseness=level)
            peak_data["levels"][str(level)] = [normalize_result(r) for r in results]

        fixture["peaks"].append(peak_data)

    with open(output_path, "w") as f:
        json.dump(fixture, f, indent=2)

    print(f"Generated golden file: {output_path}")
