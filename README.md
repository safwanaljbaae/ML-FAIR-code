# ML-FAIR: Machine-Learning-Assisted FAIR Analysis of Mean Motion Resonances

A Python framework for identifying and analyzing mean-motion resonances (MMRs) in the Solar System using N-body simulations and clustering techniques.

## Overview

ML-FAIR applies machine learning methods (KDE peak counting + OPTICS clustering) to detect and characterize mean-motion resonances from N-body integrations. The "FAIR" angle analysis helps identify resonant behavior in asteroid orbital elements.

The project consists of three independent pipelines:

| Pipeline | Description |
|----------|-------------|
| **ML-FAIR-INITIAL_CONDITIONS** | Prepare asteroid data from SBDB and generate orbital elements |
| **ML-FAIR-INTERNAL_MMR** | Analyze internal resonances (asteroids inside planet's orbit) |
| **ML-FAIR-EXTERNAL_MMR** | Analyze external resonances (asteroids outside planet's orbit) |

## Project Structure

```
ML-FAIR/
├── ml_fair/                       # Shared library (common code)
│   ├── __init__.py
│   ├── io.py                      # Input/output utilities
│   ├── analysis.py                # Core analysis functions (KDE, OPTICS)
│   └── planets.py                 # Planet data and configuration
├── ML-FAIR-INITIAL_CONDITIONS/    # Data preparation tools
│   ├── select_ast_earth_cli.py    # Filter SBDB data & query Horizons
│   └── split_rebound_cli.py       # Split data for parallel processing
├── ML-FAIR-INTERNAL_MMR/          # Internal resonance analysis
│   ├── Internal_MMR_ML_FAIR.py
│   └── ...
└── ML-FAIR-EXTERNAL_MMR/           # External resonance analysis
    ├── External_MMR_ML_FAIR.py
    └── ...
```

---

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/safwanaljbaae/ML-FAIR-code.git
cd ML-FAIR-code/
```

### 2. Create a Virtual Environment

```bash
# Create a new virtual environment
python3 -m venv .venv

# Activate it (Linux/macOS)
source .venv/bin/activate

# Activate it (Windows)
.venv\Scripts\activate
```

### 3. Install the Package

```bash
# Install ml_fair package in editable mode (includes all dependencies)
pip install -e .
```

This will install:
- **rebound** - N-body simulations
- **astropy, astroquery** - Astronomy utilities & JPL Horizons queries
- **numpy, scipy, pandas** - Scientific computing
- **matplotlib** - Visualization
- **scikit-learn** - Clustering (OPTICS)
- **tqdm** - Progress bars
- **angles** - Angle normalization

### 4. Verify Installation

```bash
python -c "from ml_fair import io, analysis, planets; print('ML-FAIR installed successfully\!')"
```

---

## Quick Start

### Prepare Initial Conditions (Optional)

If you need to prepare asteroid data from scratch:

```bash
cd ML-FAIR-INITIAL_CONDITIONS

# Download asteroid data from SBDB (NASA JPL Small-Body Database)
# Then filter asteroids by semi-major axis and query JPL Horizons
python select_ast_earth_cli.py --input-csv sbdb_query_results.csv --a-max 0.98 --sort-by a --designation-column pdes --ast-columns a,pdes --ast-output ast_VRI.txt --particles-output Particles.el --epoch 2460800.5 --location 500@10

# Split large datasets into batch directories for parallel processing
python split_rebound_cli.py --ast-file ast_VRI.txt --particles-file Particles.el --template-script Internal_MMR_ML_FAIR.py --lines-per-batch 100 --overwrite
```

### Run Internal MMR Analysis

Analyze resonances for asteroids **inside** a planet's orbit:

```bash
cd ML-FAIR-INTERNAL_MMR
python Internal_MMR_ML_FAIR.py --planet earth
```

Valid planets: `mercury`, `venus`, `earth`, `mars`, `jupiter`, `saturn`, `uranus`, `neptune`

### Run External MMR Analysis

Analyze resonances for asteroids **outside** a planet's orbit:

```bash
cd ML-FAIR-EXTERNAL_MMR
python External_MMR_ML_FAIR.py --planet mercury
```

---

## What the Code Does

For each asteroid, the pipeline:

1. **Loads data** - Reads asteroid IDs and orbital elements from `ast_VRI.txt` and `Particles.el`
2. **Transforms coordinates** - Converts orbital elements from ecliptic to invariable plane
3. **Fetches planet data** - Queries JPL Horizons for planetary orbital elements
4. **Runs N-body simulation** - Integrates Solar System + asteroid using REBOUND
5. **Computes FAIR angle** - Calculates the resonant angle σ in the invariable plane
6. **Analyzes distributions** - Uses KDE peak counting and OPTICS clustering
7. **Detects resonances** - Maps FAIR counts to resonance integers (p, q)
8. **Generates outputs** - Creates diagnostic plots and CSV summary files

### Resonant Angle Formulas

**External MMR:**
```
σ = (Ω + ω + M)_asteroid - (Ω + ω + M)_planet
a_res = a_planet × ((p+q)/p)^(2/3)
```

**Internal MMR:**
```
σ = (Ω + ω + M)_planet - (Ω + ω + M)_asteroid
a_res = a_planet × (p/(p+q))^(2/3)
```

---

## Input File Formats

### ast_VRI.txt
```
<a_input_AU> <asteroid_identifier>
```
Example:
```
0.8259 337248
0.8261 "2020 NC"
0.8269 "2004 FU162"
```

### Particles.el
```
a e i Omega omega M
```
Angles in degrees, semi-major axis in AU.

### Optional: real (ground-truth)
```
<asteroid_name> <n_clusters_sigma> <n_clusters_M>
```
Used for evaluation metrics. If absent, accuracy metrics are skipped.

---

## Output Files

For each asteroid:
- `res_forgacs_<planet>_<asteroid>.png` - Diagnostic plots
- `res_forgacs_<planet>_<asteroid>_FLAGGED.png` - Flagged case plots

Summary CSV (for resonant candidates with Δa < 0.4 AU):
- `resonant_summary_<planet>.csv` (internal)
- `resonant_summary_<planet>.csv` (external)

---

## Manual Intervention

Asteroids are **FLAGGED** when KDE peak count and OPTICS cluster count disagree. When this happens:
1. Diagnostic figures are saved
2. Script pauses and prompts for manual `nsigma` and `nM` values

This is a **semi-automatic** workflow. For batch processing, consider:
- Pre-prepared review tables
- CSV override files for flagged asteroids
- Skipping flagged objects until manual review

---

## Shared Library (ml_fair)

The `ml_fair` package provides common utilities:

```python
from ml_fair import io, analysis, planets

# Load asteroid data
ast_id, ast_a_input = io.load_asteroid_list("ast_VRI.txt")
elements = io.ler_elementos_particles()

# Get planet data
planets.get_planets()
cfg = planets.get_simulation_config("earth")

# Analyze resonances
n_peaks, peak_pos, grid, kde = analysis.count_peaks_circular(data)
a_res = analysis.resonance_a_outer(a_planet, p, pq)
```

---

## Current Limitations

- Plotting backend hard-coded as `TkAgg` (may fail on headless servers)
- Fixed filenames in working directory
- Asteroid matching is positional (by index, not name)
- JPL Horizons query requires internet access
- Multiprocessing hard-coded to `max_workers=6`

---

## License

Each subdirectory contains its own LICENSE file.

---

## Citation

Please cite the appropriate publications when using this code for academic work.
