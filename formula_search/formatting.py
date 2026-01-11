"""Output formatting utilities for formula search results."""

from typing import Dict, Any

from .constants import COARSENESS_LEVELS, LEVEL_NAMES


def format_hit(hit: Dict[str, Any]) -> str:
    """Format a single hit as a fixed-width string for display."""
    return (
        f"{hit['formula']:24s} {hit['mz']:10.4f}  {hit['charge']:+d}  "
        f"{hit['adduct']:10s} {hit['neutral_mass']:10.4f}  {hit['ppm_error']:+6.2f} ppm"
    )


def print_header() -> None:
    """Print the header row for results table."""
    print(
        f"{'Formula':24s} {'m/z':>10s}  {'z':>2s}  {'Adduct':10s} {'Neutral':>10s}  {'Error':>10s}"
    )
    print("-" * 78)


def print_results(hits: list, show_header: bool = True) -> None:
    """Print a list of hits in table format."""
    if show_header:
        print_header()
    for h in hits:
        print(format_hit(h))


def print_level_header(level: int) -> None:
    """Print a header for a coarseness level section."""
    params = COARSENESS_LEVELS[level]
    print(
        f"\n[Level {level}: {LEVEL_NAMES[level]}] "
        f"h_max={params['h_max']}, c_max={params['c_max']}, additional_o={params['additional_o']}"
    )
