"""Command-line interface for formula searching."""

import argparse

from .search import search_mz_negative
from .formatting import format_hit, print_header, print_level_header


def run_scan_all(peak_mz: float, ppm: float) -> None:
    """Run search at all coarseness levels, showing only new hits at each level."""
    prev_formulas = set()
    for level in [1, 2, 3]:
        hits = search_mz_negative(peak_mz, ppm=ppm, coarseness=level)
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
    args = parser.parse_args()

    if args.scan_all:
        run_scan_all(args.peak_mz, args.ppm)
    else:
        hits = search_mz_negative(
            args.peak_mz, ppm=args.ppm, coarseness=args.coarseness
        )
        if hits:
            print_header()
            for h in hits:
                print(format_hit(h))
        else:
            print("No matches found.")


if __name__ == "__main__":
    main()
