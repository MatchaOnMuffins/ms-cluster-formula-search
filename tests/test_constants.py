"""Tests for formula_search.constants module."""

import pytest
from formula_search.constants import (
    MASS,
    PROTON,
    CL35,
    NA23,
    K39,
    NH4,
    ADDUCTS_NEG,
    ADDUCTS_POS,
    COARSENESS_LEVELS,
    LEVEL_NAMES,
    SUPPORTED_METALS,
    get_coarseness_params,
)


class TestMassConstants:
    """Test physical mass constants."""

    def test_elemental_masses_exist(self):
        """Verify all required elemental masses are defined."""
        required = ["Y", "La", "Mn", "O", "H", "C", "tBuCOO"]
        for element in required:
            assert element in MASS, f"Missing mass for {element}"

    def test_yttrium_mass(self):
        """Yttrium monoisotopic mass should be ~88.906 Da."""
        assert abs(MASS["Y"] - 88.90584) < 1e-5

    def test_lanthanum_mass(self):
        """Lanthanum monoisotopic mass should be ~138.905 Da."""
        assert abs(MASS["La"] - 138.90547) < 1e-5

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

    def test_sodium23_mass(self):
        """Sodium-23 mass should be ~22.990 Da."""
        assert abs(NA23 - 22.98976928) < 1e-8

    def test_potassium39_mass(self):
        """Potassium-39 mass should be ~38.964 Da."""
        assert abs(K39 - 38.96370649) < 1e-8

    def test_ammonium_mass(self):
        """Ammonium (NH4) mass should be ~18.034 Da."""
        expected = 14.00307400443 + 4 * 1.00782503223  # N + 4H
        assert abs(NH4 - expected) < 1e-10
        assert 18.0 < NH4 < 18.1


class TestAdductsNegative:
    """Test negative adduct definitions."""

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


class TestAdductsPositive:
    """Test positive adduct definitions."""

    def test_positive_adducts_exist(self):
        """Verify all expected positive adducts are defined."""
        expected = ["[M+H]+", "[M]+•", "[M+Na]+", "[M+K]+", "[M+NH4]+"]
        for adduct in expected:
            assert adduct in ADDUCTS_POS, f"Missing adduct: {adduct}"

    def test_protonated_adduct(self):
        """[M+H]+ should subtract proton mass from observed."""
        assert ADDUCTS_POS["[M+H]+"] == PROTON

    def test_radical_cation_adduct(self):
        """[M]+• should have no mass adjustment."""
        assert ADDUCTS_POS["[M]+•"] == 0.0

    def test_sodiated_adduct(self):
        """[M+Na]+ should use Na-23 mass."""
        assert ADDUCTS_POS["[M+Na]+"] == NA23

    def test_potassiated_adduct(self):
        """[M+K]+ should use K-39 mass."""
        assert ADDUCTS_POS["[M+K]+"] == K39

    def test_ammoniated_adduct(self):
        """[M+NH4]+ should use NH4 mass."""
        assert ADDUCTS_POS["[M+NH4]+"] == NH4


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


class TestSupportedMetals:
    """Test supported metals configuration."""

    def test_supported_metals_contains_yttrium(self):
        """Y (yttrium) should be in SUPPORTED_METALS."""
        assert "Y" in SUPPORTED_METALS

    def test_supported_metals_contains_lanthanum(self):
        """La (lanthanum) should be in SUPPORTED_METALS."""
        assert "La" in SUPPORTED_METALS

    def test_supported_metals_is_tuple(self):
        """SUPPORTED_METALS should be a tuple."""
        assert isinstance(SUPPORTED_METALS, tuple)

    def test_all_supported_metals_have_masses(self):
        """All supported metals should have masses defined in MASS."""
        for metal in SUPPORTED_METALS:
            assert metal in MASS, f"Missing mass for supported metal {metal}"
