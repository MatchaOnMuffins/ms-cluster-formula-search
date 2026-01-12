"""Command-line interface for formula searching."""

import argparse

from .search import search_mz_negative, search_mz_positive
from .formatting import format_hit, print_header, print_level_header


def run_scan_all(peak_mz: float, ppm: float, mode: str, metal: str = "Y") -> None:
    """Run search at all coarseness levels, showing only new hits at each level."""
    search_fn = search_mz_positive if mode == "positive" else search_mz_negative
    prev_formulas = set()
    for level in [1, 2, 3]:
        hits = search_fn(peak_mz, ppm=ppm, coarseness=level, metal=metal)
        print_level_header(level)
        new_hits = [h for h in hits if h["formula"] not in prev_formulas]
        if not new_hits:
            print("  (no new hits)")
            continue
        print_header()
        for h in new_hits:
            print(format_hit(h))
            prev_formulas.add(h["formula"])


def main() -> None:
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description="Search for formula matches for a given m/z peak"
    )
    parser.add_argument("peak_mz", type=float, help="m/z value to search")
    parser.add_argument(
        "-c",
        "--coarseness",
        type=int,
        choices=[1, 2, 3],
        default=2,
        help="1=strict, 2=moderate (default), 3=loose",
    )
    parser.add_argument(
        "--scan-all", action="store_true", help="Run all coarseness levels"
    )
    parser.add_argument(
        "--ppm", type=float, default=10, help="PPM tolerance (default: 10)"
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=["positive", "negative"],
        default="negative",
        help="Ion mode: positive or negative (default: negative)",
    )
    parser.add_argument(
        "--metal",
        type=str,
        choices=["Y", "La"],
        default="Y",
        help="Metal base: Y (yttrium) or La (lanthanum) (default: Y)",
    )
    args = parser.parse_args()

    search_fn = search_mz_positive if args.mode == "positive" else search_mz_negative

    if args.scan_all:
        run_scan_all(args.peak_mz, args.ppm, args.mode, args.metal)
    else:
        hits = search_fn(args.peak_mz, ppm=args.ppm, coarseness=args.coarseness, metal=args.metal)
        if hits:
            print_header()
            for h in hits:
                print(format_hit(h))
        else:
            print("No matches found.")


if __name__ == "__main__":
    main()
