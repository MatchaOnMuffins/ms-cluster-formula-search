"""Physical constants and mass definitions for formula searching."""

import numpy as np

# Monoisotopic masses (in Daltons)
MASS = {
    "Y": 88.90584,
    "Mn": 54.938044,
    "O": 15.99491461957,
    "H": 1.00782503223,
    "C": 12.0,
}

# Derived masses
MASS["tBuCOO"] = 5 * MASS["C"] + 9 * MASS["H"] + 2 * MASS["O"]

# Particle masses
PROTON = 1.00727646688
CL35 = 34.968852682

# Adduct definitions for negative ion mode
ADDUCTS_NEG = {
    "[M−H]−": PROTON,
    "[M]−•": 0.0,
    "[M+Cl]−": CL35,
}

# Coarseness level parameters
# Controls how much extra H/C/O to allow in formula enumeration
COARSENESS_LEVELS = {
    1: {"h_max": 0, "c_max": 0, "additional_o": 0},  # strict
    2: {"h_max": 4, "c_max": 2, "additional_o": 2},  # moderate
    3: {"h_max": 10, "c_max": 5, "additional_o": 5},  # loose
}

LEVEL_NAMES = {1: "strict", 2: "moderate", 3: "loose"}


def get_coarseness_params(level: int) -> dict:
    """Get coarseness parameters for a given level.

    Args:
        level: Coarseness level (1=strict, 2=moderate, 3=loose)

    Returns:
        Dictionary with h_max, c_max, and additional_o values

    Raises:
        ValueError: If level is not 1, 2, or 3
    """
    if level not in COARSENESS_LEVELS:
        raise ValueError(f"Coarseness level must be 1, 2, or 3, got {level}")
    return COARSENESS_LEVELS[level]
