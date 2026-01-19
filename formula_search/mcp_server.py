#!/usr/bin/env python3
"""
MCP Server for Mass Spectrometry Formula Search.

This server provides tools to search for molecular formula matches
for mass spectrometry m/z peaks by enumerating possible compound combinations.
"""

import json
from typing import List, Dict, Any
from enum import Enum

from pydantic import BaseModel, Field, field_validator, ConfigDict
from mcp.server.fastmcp import FastMCP

from .search import search_mz_negative, search_mz_positive
from .constants import (
    MASS,
    ADDUCTS_NEG,
    ADDUCTS_POS,
    COARSENESS_LEVELS,
    LEVEL_NAMES,
    SUPPORTED_METALS,
)

# Initialize the MCP server
mcp = FastMCP("formula_search_mcp")


# Enums
class IonMode(str, Enum):
    """Ion mode for mass spectrometry search."""

    NEGATIVE = "negative"
    POSITIVE = "positive"


class Metal(str, Enum):
    """Supported metal bases for complexes."""

    Y = "Y"
    LA = "La"


class ResponseFormat(str, Enum):
    """Output format for tool responses."""

    MARKDOWN = "markdown"
    JSON = "json"


# Pydantic Models for Input Validation
class FormulaSearchInput(BaseModel):
    """Input model for formula search operations."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True)

    mz: float = Field(
        ..., description="Observed m/z value to search for (e.g., 1234.5)", gt=0
    )
    ppm: float = Field(
        default=10.0,
        description="Parts per million tolerance for mass matching (default: 10)",
        gt=0,
        le=100,
    )
    mode: IonMode = Field(
        default=IonMode.NEGATIVE,
        description="Ion mode: 'negative' (default) or 'positive'",
    )
    coarseness: int = Field(
        default=2,
        description="Search strictness: 1=strict, 2=moderate (default), 3=loose. Higher values allow more H/C/F/N/O in formulas.",
        ge=1,
        le=3,
    )
    metal: Metal = Field(
        default=Metal.Y,
        description="Metal base: 'Y' for yttrium (default) or 'La' for lanthanum",
    )
    max_results: int = Field(
        default=30,
        description="Maximum number of results to return per adduct/charge combination",
        ge=1,
        le=100,
    )
    response_format: ResponseFormat = Field(
        default=ResponseFormat.MARKDOWN,
        description="Output format: 'markdown' for human-readable or 'json' for machine-readable",
    )

    @field_validator("mz")
    @classmethod
    def validate_mz(cls, v: float) -> float:
        if v <= 0:
            raise ValueError("m/z value must be positive")
        return v


class ScanAllLevelsInput(BaseModel):
    """Input model for scanning all coarseness levels."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True)

    mz: float = Field(
        ..., description="Observed m/z value to search for (e.g., 1234.5)", gt=0
    )
    ppm: float = Field(
        default=10.0,
        description="Parts per million tolerance for mass matching (default: 10)",
        gt=0,
        le=100,
    )
    mode: IonMode = Field(
        default=IonMode.NEGATIVE,
        description="Ion mode: 'negative' (default) or 'positive'",
    )
    metal: Metal = Field(
        default=Metal.Y,
        description="Metal base: 'Y' for yttrium (default) or 'La' for lanthanum",
    )
    response_format: ResponseFormat = Field(
        default=ResponseFormat.MARKDOWN,
        description="Output format: 'markdown' for human-readable or 'json' for machine-readable",
    )


class GetInfoInput(BaseModel):
    """Input model for getting system information."""

    model_config = ConfigDict(str_strip_whitespace=True, validate_assignment=True)

    response_format: ResponseFormat = Field(
        default=ResponseFormat.MARKDOWN,
        description="Output format: 'markdown' for human-readable or 'json' for machine-readable",
    )


# Helper functions
def _format_hit_markdown(hit: Dict[str, Any]) -> str:
    """Format a single hit as a markdown table row."""
    return (
        f"| {hit['formula']:30s} | {hit['mz']:10.4f} | {hit['charge']:+d} | "
        f"{hit['adduct']:10s} | {hit['neutral_mass']:10.4f} | {hit['ppm_error']:+6.2f} |"
    )


def _format_results_markdown(hits: List[Dict[str, Any]], title: str = "") -> str:
    """Format search results as markdown."""
    if not hits:
        return "No matches found."

    lines = []
    if title:
        lines.append(f"## {title}")
        lines.append("")

    lines.append(f"Found {len(hits)} matching formulas:")
    lines.append("")
    lines.append(
        "| Formula | m/z | z | Adduct | Neutral Mass | PPM Error |"
    )
    lines.append(
        "|---------|-----|---|--------|--------------|-----------|"
    )

    for hit in hits:
        lines.append(_format_hit_markdown(hit))

    return "\n".join(lines)


def _search_formulas(
    mz: float,
    ppm: float,
    mode: str,
    coarseness: int,
    metal: str,
    max_hits: int,
) -> List[Dict[str, Any]]:
    """Perform formula search using the appropriate mode."""
    search_fn = search_mz_positive if mode == "positive" else search_mz_negative
    return search_fn(
        mz=mz,
        ppm=ppm,
        coarseness=coarseness,
        metal=metal,
        max_hits=max_hits,
    )


# Tool definitions
@mcp.tool(
    name="formula_search_mz",
    annotations={
        "title": "Search Formula by m/z",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def formula_search_mz(params: FormulaSearchInput) -> str:
    """
    Search for molecular formula matches for a given m/z value.

    This tool enumerates possible molecular formulas of the form:
    M_a Mn_b (tBuCOO)_c O_d H_e C_f F_g N_h

    Where M is the selected metal base (Y=yttrium or La=lanthanum),
    Mn is manganese, and tBuCOO is tert-butyl carboxylate.

    The search applies chemical constraints:
    - At least one metal (M or Mn)
    - At least one ligand (tBuCOO) or oxygen
    - Charge balance: 2*tBuCOO + O >= metals

    Args:
        params (FormulaSearchInput): Validated input parameters containing:
            - mz (float): Observed m/z value to search (e.g., 1234.5)
            - ppm (float): PPM tolerance (default: 10)
            - mode (str): Ion mode - 'negative' or 'positive' (default: negative)
            - coarseness (int): 1=strict, 2=moderate, 3=loose (default: 2)
            - metal (str): Metal base - 'Y' or 'La' (default: Y)
            - max_results (int): Maximum results per adduct (default: 30)
            - response_format (str): 'markdown' or 'json' (default: markdown)

    Returns:
        str: Search results containing matching formulas with:
            - formula: Molecular formula string
            - mz: Original m/z query value
            - charge: Charge state
            - adduct: Adduct type (e.g., [M-H]-, [M+Na]+)
            - neutral_mass: Calculated neutral mass
            - ppm_error: Parts per million error from target
    """
    hits = _search_formulas(
        mz=params.mz,
        ppm=params.ppm,
        mode=params.mode.value,
        coarseness=params.coarseness,
        metal=params.metal.value,
        max_hits=params.max_results,
    )

    if params.response_format == ResponseFormat.JSON:
        return json.dumps(
            {
                "query": {
                    "mz": params.mz,
                    "ppm": params.ppm,
                    "mode": params.mode.value,
                    "coarseness": params.coarseness,
                    "metal": params.metal.value,
                },
                "total_hits": len(hits),
                "results": hits,
            },
            indent=2,
            default=float,
        )

    title = f"Formula Search: m/z = {params.mz} ({params.mode.value} mode, {LEVEL_NAMES[params.coarseness]} coarseness)"
    return _format_results_markdown(hits, title)


@mcp.tool(
    name="formula_scan_all_levels",
    annotations={
        "title": "Scan All Coarseness Levels",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def formula_scan_all_levels(params: ScanAllLevelsInput) -> str:
    """
    Scan for formula matches at all coarseness levels (strict, moderate, loose).

    This tool runs the formula search at all three coarseness levels and shows
    only NEW hits at each level (hits not already found at stricter levels).
    This helps identify which formulas require looser search parameters.

    Args:
        params (ScanAllLevelsInput): Validated input parameters containing:
            - mz (float): Observed m/z value to search (e.g., 1234.5)
            - ppm (float): PPM tolerance (default: 10)
            - mode (str): Ion mode - 'negative' or 'positive' (default: negative)
            - metal (str): Metal base - 'Y' or 'La' (default: Y)
            - response_format (str): 'markdown' or 'json' (default: markdown)

    Returns:
        str: Search results organized by coarseness level, showing only new
             hits at each level that weren't found at stricter levels.
    """
    search_fn = (
        search_mz_positive if params.mode == IonMode.POSITIVE else search_mz_negative
    )

    all_results = {}
    prev_formulas = set()

    for level in [1, 2, 3]:
        hits = search_fn(
            mz=params.mz, ppm=params.ppm, coarseness=level, metal=params.metal.value
        )
        new_hits = [h for h in hits if h["formula"] not in prev_formulas]
        all_results[level] = {
            "level_name": LEVEL_NAMES[level],
            "params": COARSENESS_LEVELS[level],
            "new_hits": new_hits,
            "total_at_level": len(hits),
        }
        for h in new_hits:
            prev_formulas.add(h["formula"])

    if params.response_format == ResponseFormat.JSON:
        return json.dumps(
            {
                "query": {
                    "mz": params.mz,
                    "ppm": params.ppm,
                    "mode": params.mode.value,
                    "metal": params.metal.value,
                },
                "levels": {
                    str(level): {
                        "level_name": data["level_name"],
                        "params": data["params"],
                        "new_hits_count": len(data["new_hits"]),
                        "total_hits_at_level": data["total_at_level"],
                        "new_hits": data["new_hits"],
                    }
                    for level, data in all_results.items()
                },
            },
            indent=2,
            default=float,
        )

    # Markdown format
    lines = [
        f"# Formula Scan: m/z = {params.mz}",
        f"Mode: {params.mode.value}, Metal: {params.metal.value}, PPM: {params.ppm}",
        "",
    ]

    for level, data in all_results.items():
        params_str = ", ".join(f"{k}={v}" for k, v in data["params"].items())
        lines.append(f"## Level {level}: {data['level_name']} ({params_str})")
        lines.append("")

        if not data["new_hits"]:
            lines.append("*(no new hits at this level)*")
        else:
            lines.append(f"New hits: {len(data['new_hits'])}")
            lines.append("")
            lines.append(
                "| Formula | m/z | z | Adduct | Neutral Mass | PPM Error |"
            )
            lines.append(
                "|---------|-----|---|--------|--------------|-----------|"
            )
            for hit in data["new_hits"]:
                lines.append(_format_hit_markdown(hit))

        lines.append("")

    return "\n".join(lines)


@mcp.tool(
    name="formula_get_info",
    annotations={
        "title": "Get System Information",
        "readOnlyHint": True,
        "destructiveHint": False,
        "idempotentHint": True,
        "openWorldHint": False,
    },
)
async def formula_get_info(params: GetInfoInput) -> str:
    """
    Get information about the formula search system configuration.

    Returns details about supported metals, adducts, coarseness levels,
    and element masses used in the formula enumeration.

    Args:
        params (GetInfoInput): Validated input parameters containing:
            - response_format (str): 'markdown' or 'json' (default: markdown)

    Returns:
        str: System information including:
            - Supported metals (Y, La)
            - Available adducts for positive and negative modes
            - Coarseness level parameters
            - Element monoisotopic masses
    """
    info = {
        "supported_metals": list(SUPPORTED_METALS),
        "adducts": {
            "negative": {name: float(mass) for name, mass in ADDUCTS_NEG.items()},
            "positive": {name: float(mass) for name, mass in ADDUCTS_POS.items()},
        },
        "coarseness_levels": {
            str(level): {"name": LEVEL_NAMES[level], **params_dict}
            for level, params_dict in COARSENESS_LEVELS.items()
        },
        "element_masses": {elem: float(mass) for elem, mass in MASS.items()},
    }

    if params.response_format == ResponseFormat.JSON:
        return json.dumps(info, indent=2)

    # Markdown format
    lines = [
        "# Formula Search System Information",
        "",
        "## Supported Metals",
        f"- **Y** (Yttrium): {MASS['Y']:.5f} Da",
        f"- **La** (Lanthanum): {MASS['La']:.5f} Da",
        "",
        "## Adducts",
        "",
        "### Negative Ion Mode",
    ]

    for name, mass in ADDUCTS_NEG.items():
        lines.append(f"- {name}: {float(mass):.5f} Da")

    lines.extend(["", "### Positive Ion Mode"])

    for name, mass in ADDUCTS_POS.items():
        lines.append(f"- {name}: {float(mass):.5f} Da")

    lines.extend(["", "## Coarseness Levels", ""])

    for level, params_dict in COARSENESS_LEVELS.items():
        params_str = ", ".join(f"{k}={v}" for k, v in params_dict.items())
        lines.append(f"- **Level {level} ({LEVEL_NAMES[level]})**: {params_str}")

    lines.extend(["", "## Element Masses (Monoisotopic)", ""])

    for elem, mass in MASS.items():
        if elem != "tBuCOO":
            lines.append(f"- **{elem}**: {float(mass):.5f} Da")
        else:
            lines.append(f"- **{elem}** (5C + 9H + 2O): {float(mass):.5f} Da")

    return "\n".join(lines)


if __name__ == "__main__":
    mcp.run()
