import numpy as np

# Monoisotopic masses
MASS = {
    "Y": 88.90584,
    "Mn": 54.938044,
    "O": 15.99491461957,
    "H": 1.00782503223,
    "C": 12.0,
}

MASS["tBuCOO"] = 5 * MASS["C"] + 9 * MASS["H"] + 2 * MASS["O"]

PROTON = 1.00727646688
CL35 = 34.968852682

ADDUCTS_NEG = {
    "[M−H]−": PROTON,
    "[M]−•": 0.0,
    "[M+Cl]−": CL35,
}


def within_ppm(m, target, ppm):
    return abs(m - target) <= target * ppm * 1e-6


COARSENESS_LEVELS = {
    1: {"h_max": 0, "c_max": 0, "additional_o": 0},
    2: {"h_max": 4, "c_max": 2, "additional_o": 2},
    3: {"h_max": 10, "c_max": 5, "additional_o": 5},
}


def get_coarseness_params(level):
    if level not in COARSENESS_LEVELS:
        raise ValueError(f"Coarseness level must be 1, 2, or 3, got {level}")
    return COARSENESS_LEVELS[level]


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
    Enumerate formulas Y_a Mn_b (tBuCOO)_c O_d H_e C_f matching target_mass within ppm.

    Constraints:
      - At least one metal (Y or Mn)
      - At least one ligand (tBuCOO) or oxygen
      - If tBuCOO present: 2*tBuCOO + O >= metals (charge balance)

    Coarseness (1-3) controls extra H/C/O: 1=strict, 2=moderate, 3=loose.
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

    valid = ~((y == 0) & (mn == 0))
    valid &= ~((k == 0) & (o == 0))
    valid &= (k == 0) | ((2 * k + o) >= (mn + y))

    y, mn, k, o, h, c = y[valid], mn[valid], k[valid], o[valid], h[valid], c[valid]

    masses = (
        y * MASS["Y"]
        + mn * MASS["Mn"]
        + k * MASS["tBuCOO"]
        + o * MASS["O"]
        + h * MASS["H"]
        + c * MASS["C"]
    )

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


def search_mz_negative(
    mz,
    ppm=5,
    charges=(-1,),
    adducts=ADDUCTS_NEG,
    y_max=2,
    mn_max=5,
    tbu_max=11,
    o_max=5,
    h_max=None,
    c_max=None,
    coarseness=2,
    max_hits=30,
):
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


LEVEL_NAMES = {1: "strict", 2: "moderate", 3: "loose"}


def format_hit(hit):
    return (
        f"{hit['formula']:24s} {hit['mz']:10.4f}  {hit['charge']:+d}  "
        f"{hit['adduct']:10s} {hit['neutral_mass']:10.4f}  {hit['ppm_error']:+6.2f} ppm"
    )


def print_header():
    print(f"{'Formula':24s} {'m/z':>10s}  {'z':>2s}  {'Adduct':10s} {'Neutral':>10s}  {'Error':>10s}")
    print("-" * 78)


def run_scan_all(peak_mz, ppm):
    prev_formulas = set()
    for level in [1, 2, 3]:
        hits = search_mz_negative(peak_mz, ppm=ppm, coarseness=level)
        params = COARSENESS_LEVELS[level]
        print(f"\n[Level {level}: {LEVEL_NAMES[level]}] "
              f"h_max={params['h_max']}, c_max={params['c_max']}, additional_o={params['additional_o']}")
        new_hits = [h for h in hits if h["formula"] not in prev_formulas]
        if not new_hits:
            print("  (no new hits)")
            continue
        print_header()
        for h in new_hits:
            print(format_hit(h))
            prev_formulas.add(h["formula"])


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Search for formula matches for a given m/z peak")
    parser.add_argument("peak_mz", type=float, help="m/z value to search")
    parser.add_argument("-c", "--coarseness", type=int, choices=[1, 2, 3], default=2,
                        help="1=strict, 2=moderate (default), 3=loose")
    parser.add_argument("--scan-all", action="store_true", help="Run all coarseness levels")
    parser.add_argument("--ppm", type=float, default=10, help="PPM tolerance (default: 10)")
    args = parser.parse_args()

    if args.scan_all:
        run_scan_all(args.peak_mz, args.ppm)
    else:
        hits = search_mz_negative(args.peak_mz, ppm=args.ppm, coarseness=args.coarseness)
        if hits:
            print_header()
            for h in hits:
                print(format_hit(h))
        else:
            print("No matches found.")


if __name__ == "__main__":
    main()
