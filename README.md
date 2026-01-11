# ms-cluster-formula-search

## Overview

This project enumerates molecular formulas of the form

$$
\mathrm{Y_{a}Mn_{b}(tBuCOO)_{c}O_{d}H_{e}C_{f}}
$$

that correspond to observed mass spectrometry peaks, matching a specified target mass within a given parts-per-million (ppm) tolerance.

## Quick Start

```bash
# Search for a specific m/z value (e.g., 1234.5)
./search.py 1234.5

# Adjust PPM tolerance (default: 10)
./search.py 1234.5 --ppm 5

# Change coarseness level: 1=strict, 2=moderate (default), 3=loose
./search.py 1234.5 -c 1

# Scan all coarseness levels at once
./search.py 1234.5 --scan-all
```
