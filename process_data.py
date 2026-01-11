import numpy as np

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
# Coarseness levels for secondary adducts (C, H)
# ----------------------------
COARSENESS_LEVELS = {
    1: {"h_max": 0, "c_max": 0, "additional_o": 0},   # strict: no extra C, H
    2: {"h_max": 4, "c_max": 2, "additional_o": 2},   # moderate: a few extra C, H
    3: {"h_max": 10, "c_max": 5, "additional_o": 5}, # loose: up to 10 extra H, 5 extra C
}


def get_coarseness_params(level):
    """Get h_max and c_max for a given coarseness level (1-3)."""
    if level not in COARSENESS_LEVELS:
        raise ValueError(f"Coarseness level must be 1, 2, or 3, got {level}")
    return COARSENESS_LEVELS[level]


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
    h_max=None,
    c_max=None,
    coarseness=2,
):
    """
    Enumerate neutral formulas:
        Y_x Mn_y (tBuCOO)_k O_m H_n C_c
    matching target_mass within ppm.

    coarseness: 1-3, controls h_max, c_max, and additional_o if not explicitly set
        1 = strict (no extra C, H, O)
        2 = moderate (a few extra C, H, O)
        3 = loose (up to 10 extra C, H and 5 extra O)
    """
    params = get_coarseness_params(coarseness)
    if h_max is None:
        h_max = params["h_max"]
    if c_max is None:
        c_max = params["c_max"]
    o_max_total = o_max + params["additional_o"]

    # Generate all index combinations using meshgrid
    y_vals = np.arange(y_max + 1)
    mn_vals = np.arange(mn_max + 1)
    k_vals = np.arange(tbu_max + 1)
    o_vals = np.arange(o_max_total + 1)
    h_vals = np.arange(h_max + 1)
    c_vals = np.arange(c_max + 1)

    y, mn, k, o, h, c = np.meshgrid(y_vals, mn_vals, k_vals, o_vals, h_vals, c_vals, indexing='ij')
    y = y.ravel()
    mn = mn.ravel()
    k = k.ravel()
    o = o.ravel()
    h = h.ravel()
    c = c.ravel()

    # Apply filters vectorized
    # Filter 1: y == 0 and mn == 0 -> exclude
    valid = ~((y == 0) & (mn == 0))
    # Filter 2: k == 0 and o == 0 -> exclude
    valid &= ~((k == 0) & (o == 0))
    # Filter 3: if k > 0, need (2*k + o) >= (mn + y)
    valid &= (k == 0) | ((2 * k + o) >= (mn + y))

    # Apply filter
    y, mn, k, o, h, c = y[valid], mn[valid], k[valid], o[valid], h[valid], c[valid]

    # Compute masses vectorized
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

    y, mn, k, o, h, c = y[ppm_mask], mn[ppm_mask], k[ppm_mask], o[ppm_mask], h[ppm_mask], c[ppm_mask]
    masses = masses[ppm_mask]
    ppm_errors = (masses - target_mass) / target_mass * 1e6

    # Sort by absolute ppm error
    sort_idx = np.argsort(np.abs(ppm_errors))

    # Build result list
    hits = []
    for i in sort_idx:
        hits.append({
            "formula": f"Y{y[i]}Mn{mn[i]}(tBuCOO){k[i]}O{o[i]}H{h[i]}C{c[i]}",
            "mass": masses[i],
            "ppm_error": ppm_errors[i],
            "counts": (int(y[i]), int(mn[i]), int(k[i]), int(o[i]), int(h[i]), int(c[i])),
        })

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
    h_max=None,
    c_max=None,
    coarseness=2,
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
                    c_max=c_max,
                    coarseness=coarseness,
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
    import argparse

    parser = argparse.ArgumentParser(description="Search for formula matches for a given m/z peak")
    parser.add_argument("peak_mz", type=float, help="The m/z value to search")
    parser.add_argument(
        "-c", "--coarseness",
        type=int,
        choices=[1, 2, 3],
        default=2,
        help="Coarseness level for C/H enumeration: 1=strict (no extra C,H), 2=moderate, 3=loose (up to 10 extra C,H)"
    )
    parser.add_argument(
        "--scan-all",
        action="store_true",
        help="Run all coarseness levels (1, 2, 3) and show results for each"
    )
    parser.add_argument("--ppm", type=float, default=10, help="PPM tolerance (default: 10)")
    args = parser.parse_args()

    peaks = [args.peak_mz]

    if args.scan_all:
        level_names = {1: "strict", 2: "moderate", 3: "loose"}
        prev_formulas = set()
        for level in [1, 2, 3]:
            hits = scan_peaks_tBuCOO_YMn_negative(
                peaks,
                ppm=args.ppm,
                charges=(-1,),
                y_max=2,
                mn_max=5,
                tbu_max=11,
                o_max=5,
                coarseness=level,
            )
            params = COARSENESS_LEVELS[level]
            print(f"\n=== Level {level} ({level_names[level]}): h_max={params['h_max']}, c_max={params['c_max']}, additional_o={params['additional_o']} ===")
            new_hits = [h for h in hits if h['formula'] not in prev_formulas]
            if not new_hits:
                print("  (no new hits)")
            for h in new_hits:
                print(f"  {h['formula']:20s} mz={h['mz']:.4f} "
                      f"charge={h['charge']} adduct={h['adduct']} "
                      f"neutral_mass={h['neutral_mass']:.4f} ppm_error={h['ppm_error']:.2f}")
                prev_formulas.add(h['formula'])
    else:
        hits = scan_peaks_tBuCOO_YMn_negative(
            peaks,
            ppm=args.ppm,
            charges=(-1,),
            y_max=2,
            mn_max=5,
            tbu_max=11,
            o_max=5,
            coarseness=args.coarseness,
        )
        for h in hits:
            print(f"{h['formula']:20s} mz={h['mz']:.4f} "
                  f"charge={h['charge']} adduct={h['adduct']} "
                  f"neutral_mass={h['neutral_mass']:.4f} ppm_error={h['ppm_error']:.2f}")
