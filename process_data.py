# ----------------------------
# Exact monoisotopic masses
# ----------------------------
MASS = {
    "Y":  88.90584,             # 89Y
    "Mn": 54.938044,
    "O":  15.99491461957,
    "H":  1.00782503223,        # neutral H atom (for neutral formula mass)
    "C":  12.0,
}

# tBuCOO = C5H9O2
MASS["tBuCOO"] = (
    5 * MASS["C"]
    + 9 * MASS["H"]
    + 2 * MASS["O"]
)

# For converting observed m/z -> neutral mass in NEGATIVE mode
PROTON = 1.00727646688          # H+ mass (used for [M-H]- relationship)
CL35   = 34.968852682           # 35Cl- mass (common chloride adduct)

ADDUCTS_NEG = {
    # core negative-mode ions
    "[M−H]−":    PROTON,     # neutral = mz + H
    "[M]−•":     0.0,
    "[M+Cl]−":   CL35
}

def within_ppm(m, target, ppm):
    return abs(m - target) <= target * ppm * 1e-6


# ----------------------------
# Building-block enumerator
# ----------------------------
def enumerate_tBuCOO_YMn(
    target_mass,
    ppm=5,
    y_max=2,
    mn_max=5,
    tbu_max=11,
    o_max=5,
    h_max=6,
):
    """
    Enumerate neutral formulas:
        Y_x Mn_y (tBuCOO)_k O_m H_n
    matching target_mass within ppm.
    """
    hits = []
    for y in range(y_max + 1):
        for mn in range(mn_max + 1):
            if y == 0 and mn == 0:
                continue

            for k in range(tbu_max + 1):
                for o in range(o_max + 1):
                    for h in range(h_max + 1):

                        # minimal sanity
                        if k == 0 and o == 0:
                            continue

                        # crude "enough O donors for metals" constraint:
                        # each tBuCOO contributes 2 O donors, plus free O.
                        if k > 0 and (2*k + o) < (mn + y):
                            continue

                        mass = (
                            y  * MASS["Y"]
                            + mn * MASS["Mn"]
                            + k  * MASS["tBuCOO"]
                            + o  * MASS["O"]
                            + h  * MASS["H"]
                        )

                        if within_ppm(mass, target_mass, ppm):
                            hits.append({
                                "formula": f"Y{y}Mn{mn}(tBuCOO){k}O{o}H{h}",
                                "mass": mass,
                                "ppm_error": (mass - target_mass) / target_mass * 1e6,
                                "counts": (y, mn, k, o, h),
                            })

    hits.sort(key=lambda x: abs(x["ppm_error"]))
    return hits


# ----------------------------
# Scan ~2k peaks (negative mode)
# ----------------------------
def scan_peaks_tBuCOO_YMn_negative(
    peaks,
    ppm=5,
    charges=(-1,),              # usually just -1 for MALDI-/ESI-; add -2 if you truly expect it
    adducts=ADDUCTS_NEG,
    y_max=2,
    mn_max=5,
    tbu_max=11,
    o_max=5,
    h_max=6,
    max_hits_per_peak=30,       # keeps output manageable
):
    """
    peaks: iterable of (mz, intensity)

    For negative mode we use:
        neutral_mass = mz * abs(charge) + adduct_mass

    Returns list of match dicts sorted by |ppm_error| then intensity.
    """
    results = []

    for mz in peaks:
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
                )

                for h in hits[:max_hits_per_peak]:
                    results.append({
                        "mz": mz,
                        "charge": z,
                        "adduct": ad_name,
                        "neutral_mass": neutral_mass,
                        **h,
                    })

    results.sort(key=lambda r: abs(r["ppm_error"]))
    return results



if __name__ == "__main__":
    # user can enter a peak in stdin, e.g. "1234.5678"
    # for quick testing
    import sys

    if len(sys.argv) != 2:
        print("Usage: python process_data.py <peak_mz>")
        sys.exit(1)

    peak_mz = float(sys.argv[1])
    peaks = [peak_mz]  # dummy intensity
    hits = scan_peaks_tBuCOO_YMn_negative(
        peaks,
        ppm=10,
        charges=(-1,),
        y_max=2,
        mn_max=5,
        tbu_max=11,
        o_max=5,
        h_max=6,
    )
    for h in hits:
        # format this, formula first
        print(f"{h['formula']:20s} mz={h['mz']:.4f} "
              f"charge={h['charge']} adduct={h['adduct']} "
              f"neutral_mass={h['neutral_mass']:.4f} ppm_error={h['ppm_error']:.2f}")