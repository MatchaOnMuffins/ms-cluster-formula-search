"""Mass spectrometry m/z search functions."""

from typing import List, Dict, Any, Tuple

from .constants import ADDUCTS_NEG
from .enumeration import enumerate_tBuCOO_YMn


def search_mz_negative(
    mz: float,
    ppm: float = 5,
    charges: Tuple[int, ...] = (-1,),
    adducts: Dict[str, float] = None,
    y_max: int = 2,
    mn_max: int = 5,
    tbu_max: int = 11,
    o_max: int = 5,
    h_max: int = None,
    c_max: int = None,
    coarseness: int = 2,
    max_hits: int = 30,
) -> List[Dict[str, Any]]:
    """
    Search for formula matches for a given m/z in negative ion mode.

    Args:
        mz: Observed m/z value
        ppm: Parts per million tolerance
        charges: Tuple of charge states to consider
        adducts: Dictionary of adduct names to mass adjustments
        y_max: Maximum yttrium count
        mn_max: Maximum manganese count
        tbu_max: Maximum tert-butyl carboxylate count
        o_max: Base maximum oxygen count
        h_max: Maximum hydrogen count (overrides coarseness)
        c_max: Maximum carbon count (overrides coarseness)
        coarseness: 1=strict, 2=moderate, 3=loose
        max_hits: Maximum hits per adduct/charge combination

    Returns:
        List of result dictionaries sorted by absolute ppm error
    """
    if adducts is None:
        adducts = ADDUCTS_NEG

    results = []
    mz = float(mz)

    for z in charges:
        for ad_name, ad_mass in adducts.items():
            neutral_mass = mz * abs(z) + ad_mass
            if neutral_mass <= 0:
                continue

            hits = enumerate_tBuCOO_YMn(
                neutral_mass,
                ppm=ppm,
                y_max=y_max,
                mn_max=mn_max,
                tbu_max=tbu_max,
                o_max=o_max,
                h_max=h_max,
                c_max=c_max,
                coarseness=coarseness,
            )

            for h in hits[:max_hits]:
                results.append(
                    {
                        "mz": mz,
                        "charge": z,
                        "adduct": ad_name,
                        "neutral_mass": neutral_mass,
                        **h,
                    }
                )

    results.sort(key=lambda r: abs(r["ppm_error"]))
    return results
