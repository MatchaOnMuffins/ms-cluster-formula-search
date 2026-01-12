"""Golden file tests - validate search results against known expected outputs.

These tests read expected outputs from JSON fixtures and verify that the
search functions produce matching results. This ensures reproducibility
and catches regressions in the formula enumeration algorithm.

Also includes integration tests that verify the complete workflow.
"""

import json
import io
import sys
import pytest
from pathlib import Path
from typing import Any

from formula_search import (
    search_mz_negative,
    search_mz_positive,
    format_hit,
    run_scan_all,
)


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
    """Normalize a result dict for comparison."""
    return {
        "formula": result["formula"],
        "mz": result["mz"],
        "charge": result["charge"],
        "adduct": result["adduct"],
        "neutral_mass": round(result["neutral_mass"], 6),
        "ppm_error": round(result["ppm_error"], 4),
        "counts": list(result["counts"]),
    }


class TestGoldenFixture:
    """Test that the golden fixture file is valid."""

    def test_fixture_file_exists(self):
        """Verify the fixture file exists and is valid JSON."""
        fixture_path = FIXTURES_DIR / "expected_outputs.json"
        assert fixture_path.exists(), f"Fixture file not found: {fixture_path}"

        with open(fixture_path) as f:
            data = json.load(f)

        assert "peaks" in data
        assert len(data["peaks"]) > 0

    def test_fixture_structure(self):
        """Verify fixture has expected structure."""
        expected_data = load_expected_outputs()
        for peak in expected_data["peaks"]:
            assert "mz" in peak
            assert "ppm" in peak
            assert "levels" in peak
            assert "mode" in peak
            assert "metal" in peak
            assert "polarity" in peak
            assert all(level in peak["levels"] for level in ["1", "2", "3"])

    def test_all_modes_present(self):
        """Verify all 4 modes have fixture data."""
        expected_data = load_expected_outputs()
        modes_in_fixture = {p["mode"] for p in expected_data["peaks"]}
        expected_modes = {"y_negative", "y_positive", "la_negative", "la_positive"}
        assert expected_modes.issubset(modes_in_fixture), (
            f"Missing modes: {expected_modes - modes_in_fixture}"
        )


class TestGoldenOutputs:
    """Test search results against golden (expected) outputs from fixtures."""

    @pytest.fixture
    def expected_data(self) -> dict:
        """Load the expected outputs fixture."""
        return load_expected_outputs()

    @pytest.mark.parametrize("level", [1, 2, 3])
    @pytest.mark.parametrize(
        "mode,search_fn,metal",
        MODES,
        ids=["y_negative", "y_positive", "la_negative", "la_positive"],
    )
    def test_peak_1519_matches_expected(
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


class TestResultOrdering:
    """Test that results are returned in the expected order."""

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


class TestLevelProgression:
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

        results_2 = search_fn(peak["mz"], ppm=peak["ppm"], coarseness=2, metal=metal)
        results_3 = search_fn(peak["mz"], ppm=peak["ppm"], coarseness=3, metal=metal)

        assert len(results_3) >= len(results_2), (
            f"{mode}: Level 3 should find at least as many results as level 2"
        )

    @pytest.mark.parametrize(
        "mode,search_fn,metal",
        MODES,
        ids=["y_negative", "y_positive", "la_negative", "la_positive"],
    )
    def test_level2_includes_level1_formulas(
        self, expected_data, mode, search_fn, metal
    ):
        """Verify level 2 includes all level 1 formulas."""
        peak = get_peak_by_mode(expected_data, mode)

        results_1 = search_fn(peak["mz"], ppm=peak["ppm"], coarseness=1, metal=metal)
        results_2 = search_fn(peak["mz"], ppm=peak["ppm"], coarseness=2, metal=metal)

        formulas_1 = {r["formula"] for r in results_1}
        formulas_2 = {r["formula"] for r in results_2}

        assert formulas_1.issubset(formulas_2), (
            f"{mode}: Level 2 should include all level 1 formulas"
        )


class TestIntegrationWorkflow:
    """Integration tests for complete search and formatting workflow."""

    @pytest.fixture
    def example_mz(self):
        """The example m/z value from user's output."""
        return 1519.1540

    def test_full_search_and_format_workflow(self, example_mz):
        """Test complete search and formatting workflow."""
        results = search_mz_negative(example_mz, ppm=10, coarseness=2)

        assert len(results) > 0

        for r in results:
            formatted = format_hit(r)
            assert isinstance(formatted, str)
            assert len(formatted) > 40
            assert r["formula"] in formatted

    def test_scan_all_workflow(self, example_mz):
        """Test the run_scan_all function produces output."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            run_scan_all(example_mz, ppm=10, mode="negative")
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()

        # Should have level headers
        assert "Level 1" in output
        assert "Level 2" in output
        assert "Level 3" in output

        # Should have some formulas
        assert "Y" in output
        assert "Mn" in output

    def test_scan_all_positive_mode(self):
        """Test run_scan_all with positive mode."""
        captured = io.StringIO()
        sys.stdout = captured
        try:
            run_scan_all(500.0, ppm=50, mode="positive")
        finally:
            sys.stdout = sys.__stdout__

        output = captured.getvalue()
        assert "Level" in output


# Utility function for regenerating golden file
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
