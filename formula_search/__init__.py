"""Formula search package for mass spectrometry m/z peak identification."""

from .constants import (
    MASS,
    PROTON,
    CL35,
    ADDUCTS_NEG,
    COARSENESS_LEVELS,
    LEVEL_NAMES,
    get_coarseness_params,
)
from .enumeration import enumerate_tBuCOO_YMn, within_ppm
from .search import search_mz_negative
from .formatting import format_hit, print_header, print_results, print_level_header
from .cli import main, run_scan_all

__all__ = [
    # Constants
    "MASS",
    "PROTON",
    "CL35",
    "ADDUCTS_NEG",
    "COARSENESS_LEVELS",
    "LEVEL_NAMES",
    "get_coarseness_params",
    # Enumeration
    "enumerate_tBuCOO_YMn",
    "within_ppm",
    # Search
    "search_mz_negative",
    # Formatting
    "format_hit",
    "print_header",
    "print_results",
    "print_level_header",
    # CLI
    "main",
    "run_scan_all",
]
