"""Formula search package for mass spectrometry m/z peak identification."""

from .constants import (
    MASS,
    PROTON,
    CL35,
    NA23,
    K39,
    NH4,
    ADDUCTS_NEG,
    ADDUCTS_POS,
    COARSENESS_LEVELS,
    LEVEL_NAMES,
    SUPPORTED_METALS,
    get_coarseness_params,
)
from .enumeration import enumerate_tBuCOO_YMn, within_ppm
from .search import search_mz_negative, search_mz_positive
from .formatting import format_hit, print_header, print_results, print_level_header
from .cli import main, run_scan_all

__all__ = [
    # Constants
    "MASS",
    "PROTON",
    "CL35",
    "NA23",
    "K39",
    "NH4",
    "ADDUCTS_NEG",
    "ADDUCTS_POS",
    "COARSENESS_LEVELS",
    "LEVEL_NAMES",
    "SUPPORTED_METALS",
    "get_coarseness_params",
    # Enumeration
    "enumerate_tBuCOO_YMn",
    "within_ppm",
    # Search
    "search_mz_negative",
    "search_mz_positive",
    # Formatting
    "format_hit",
    "print_header",
    "print_results",
    "print_level_header",
    # CLI
    "main",
    "run_scan_all",
]
