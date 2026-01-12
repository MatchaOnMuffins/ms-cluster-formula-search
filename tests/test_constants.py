"""Tests for formula_search.constants module."""

import pytest
from formula_search.constants import (
    MASS,
    PROTON,
    CL35,
    ADDUCTS_NEG,
    COARSENESS_LEVELS,
    LEVEL_NAMES,
    get_coarseness_params,
)


class TestMassConstants:
    """Test physical mass constants."""

    def test_elemental_masses_exist(self):
        """Verify all required elemental masses are defined."""
        required = ["Y", "Mn", "O", "H", "C", "tBuCOO"]
        for element in required:
            assert element in MASS, f"Missing mass for {element}"

    def test_yttrium_mass(self):
        """Yttrium monoisotopic mass should be ~88.906 Da."""
        assert abs(MASS["Y"] - 88.90584) < 1e-5

    def test_manganese_mass(self):
        """Manganese monoisotopic mass should be ~54.938 Da."""
        assert abs(MASS["Mn"] - 54.938044) < 1e-6

    def test_oxygen_mass(self):
        """Oxygen monoisotopic mass should be ~15.995 Da."""
        assert abs(MASS["O"] - 15.99491461957) < 1e-10

    def test_hydrogen_mass(self):
        """Hydrogen monoisotopic mass should be ~1.008 Da."""
        assert abs(MASS["H"] - 1.00782503223) < 1e-10

    def test_carbon_mass(self):
        """Carbon-12 mass should be exactly 12.0 Da."""
        assert MASS["C"] == 12.0

    def test_tBuCOO_mass_calculation(self):
        """tBuCOO mass should equal 5*C + 9*H + 2*O."""
        expected = 5 * MASS["C"] + 9 * MASS["H"] + 2 * MASS["O"]
        assert abs(MASS["tBuCOO"] - expected) < 1e-10

    def test_tBuCOO_mass_value(self):
        """tBuCOO group should be ~101.06 Da."""
        # (CH3)3C-COO- = C5H9O2
        assert 101 < MASS["tBuCOO"] < 102


class TestParticleMasses:
    """Test particle masses."""

    def test_proton_mass(self):
        """Proton mass should be ~1.007 Da."""
        assert abs(PROTON - 1.00727646688) < 1e-10

    def test_chlorine35_mass(self):
        """Chlorine-35 mass should be ~34.969 Da."""
        assert abs(CL35 - 34.968852682) < 1e-9


class TestAdducts:
    """Test adduct definitions."""

    def test_negative_adducts_exist(self):
        """Verify all expected negative adducts are defined."""
        expected = ["[M−H]−", "[M]−•", "[M+Cl]−"]
        for adduct in expected:
            assert adduct in ADDUCTS_NEG, f"Missing adduct: {adduct}"

    def test_deprotonated_adduct(self):
        """[M-H]- should add proton mass to neutral."""
        assert ADDUCTS_NEG["[M−H]−"] == PROTON

    def test_radical_anion_adduct(self):
        """[M]-• should have no mass adjustment."""
        assert ADDUCTS_NEG["[M]−•"] == 0.0

    def test_chloride_adduct(self):
        """[M+Cl]- should add Cl-35 mass."""
        assert ADDUCTS_NEG["[M+Cl]−"] == CL35


class TestCoarsenessLevels:
    """Test coarseness level parameters."""

    def test_all_levels_exist(self):
        """Levels 1, 2, and 3 should all be defined."""
        assert set(COARSENESS_LEVELS.keys()) == {1, 2, 3}

    def test_level1_strict(self):
        """Level 1 (strict) should have zero additional atoms."""
        params = COARSENESS_LEVELS[1]
        assert params["h_max"] == 0
        assert params["c_max"] == 0
        assert params["additional_o"] == 0

    def test_level2_moderate(self):
        """Level 2 (moderate) should have moderate atom limits."""
        params = COARSENESS_LEVELS[2]
        assert params["h_max"] == 4
        assert params["c_max"] == 2
        assert params["additional_o"] == 2

    def test_level3_loose(self):
        """Level 3 (loose) should have higher atom limits."""
        params = COARSENESS_LEVELS[3]
        assert params["h_max"] == 10
        assert params["c_max"] == 5
        assert params["additional_o"] == 5

    def test_level_names(self):
        """Level names should be correctly mapped."""
        assert LEVEL_NAMES[1] == "strict"
        assert LEVEL_NAMES[2] == "moderate"
        assert LEVEL_NAMES[3] == "loose"


class TestGetCoarsenessParams:
    """Test get_coarseness_params function."""

    def test_valid_levels(self):
        """Valid levels should return correct parameters."""
        for level in [1, 2, 3]:
            params = get_coarseness_params(level)
            assert params == COARSENESS_LEVELS[level]

    def test_invalid_level_zero(self):
        """Level 0 should raise ValueError."""
        with pytest.raises(ValueError, match="must be 1, 2, or 3"):
            get_coarseness_params(0)

    def test_invalid_level_four(self):
        """Level 4 should raise ValueError."""
        with pytest.raises(ValueError, match="must be 1, 2, or 3"):
            get_coarseness_params(4)

    def test_invalid_level_negative(self):
        """Negative level should raise ValueError."""
        with pytest.raises(ValueError, match="must be 1, 2, or 3"):
            get_coarseness_params(-1)
