"""Formula enumeration using vectorized numpy operations."""

import numpy as np
from typing import List, Dict, Any

from .constants import MASS, get_coarseness_params


def within_ppm(m: float, target: float, ppm: float) -> bool:
    """Check if mass m is within ppm tolerance of target."""
    return abs(m - target) <= target * ppm * 1e-6


def enumerate_tBuCOO_YMn(
    target_mass: float,
    ppm: float = 5,
    y_max: int = 2,
    mn_max: int = 5,
    tbu_max: int = 11,
    o_max: int = 5,
    h_max: int = None, # pyright: ignore[reportArgumentType]
    c_max: int = None, # pyright: ignore[reportArgumentType]
    coarseness: int = 2,
) -> List[Dict[str, Any]]:
    """
    Enumerate formulas Y_a Mn_b (tBuCOO)_c O_d H_e C_f matching target_mass within ppm.

    Constraints:
      - At least one metal (Y or Mn)
      - At least one ligand (tBuCOO) or oxygen
      - If tBuCOO present: 2*tBuCOO + O >= metals (charge balance)

    Args:
        target_mass: Target neutral mass to match
        ppm: Parts per million tolerance
        y_max: Maximum yttrium count
        mn_max: Maximum manganese count
        tbu_max: Maximum tert-butyl carboxylate count
        o_max: Base maximum oxygen count
        h_max: Maximum hydrogen count (overrides coarseness)
        c_max: Maximum carbon count (overrides coarseness)
        coarseness: 1=strict, 2=moderate, 3=loose

    Returns:
        List of hit dictionaries sorted by absolute ppm error
    """
    params = get_coarseness_params(coarseness)
    if h_max is None:
        h_max = params["h_max"]
    if c_max is None:
        c_max = params["c_max"]
    o_max_total = o_max + params["additional_o"]

    y_vals = np.arange(y_max + 1)
    mn_vals = np.arange(mn_max + 1)
    k_vals = np.arange(tbu_max + 1)
    o_vals = np.arange(o_max_total + 1)
    h_vals = np.arange(h_max + 1)
    c_vals = np.arange(c_max + 1)

    y, mn, k, o, h, c = np.meshgrid(
        y_vals, mn_vals, k_vals, o_vals, h_vals, c_vals, indexing="ij"
    )
    y = y.ravel()
    mn = mn.ravel()
    k = k.ravel()
    o = o.ravel()
    h = h.ravel()
    c = c.ravel()

    # Apply constraints
    valid = ~((y == 0) & (mn == 0))  # At least one metal
    valid &= ~((k == 0) & (o == 0))  # At least one ligand or oxygen
    valid &= (k == 0) | ((2 * k + o) >= (mn + y))  # Charge balance

    y, mn, k, o, h, c = y[valid], mn[valid], k[valid], o[valid], h[valid], c[valid]

    # Calculate masses
    masses = (
        y * MASS["Y"]
        + mn * MASS["Mn"]
        + k * MASS["tBuCOO"]
        + o * MASS["O"]
        + h * MASS["H"]
        + c * MASS["C"]
    )

    # Filter by ppm tolerance
    tol = target_mass * ppm * 1e-6
    ppm_mask = np.abs(masses - target_mass) <= tol

    y, mn, k, o, h, c = (
        y[ppm_mask],
        mn[ppm_mask],
        k[ppm_mask],
        o[ppm_mask],
        h[ppm_mask],
        c[ppm_mask],
    )
    masses = masses[ppm_mask]
    ppm_errors = (masses - target_mass) / target_mass * 1e6

    # Sort by absolute ppm error
    sort_idx = np.argsort(np.abs(ppm_errors))

    hits = []
    for i in sort_idx:
        hits.append(
            {
                "formula": f"Y{y[i]}Mn{mn[i]}(tBuCOO){k[i]}O{o[i]}H{h[i]}C{c[i]}",
                "mass": masses[i],
                "ppm_error": ppm_errors[i],
                "counts": (
                    int(y[i]),
                    int(mn[i]),
                    int(k[i]),
                    int(o[i]),
                    int(h[i]),
                    int(c[i]),
                ),
            }
        )

    return hits
