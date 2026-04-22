# FAIR Initial-Condition Preparation

This branch contains a small command-line workflow for preparing initial-condition files for FAIR-based simulations from a Small-Body Database (SBDB) query table.

The workflow is split into two stages:

1. **Filter and query orbital elements** with `select_ast_earth_cli.py`.
2. **Split the resulting files into batch directories** with `split_rebound_cli.py`.

These scripts are intended to support a repository branch focused on generating simulation-ready initial conditions before downstream FAIR / REBOUND-style integrations are launched.

## Included files

- `select_ast_earth_cli.py` — reads an SBDB-style CSV table, filters objects by semi-major axis, queries JPL Horizons, and writes simulation input files.
- `split_rebound_cli.py` — partitions the generated files into numbered batch directories and copies a simulation template into each directory.
- `run_select_earth` — minimal wrapper showing a typical call to the selection script.
- `run_split_rebound_cli` — minimal wrapper showing a typical call to the batch-splitting script.

## What the workflow produces

After a successful run of `select_ast_earth_cli.py`, the main outputs are:

- `ast_VRI.txt` — filtered asteroid list
- `Particles.el` — orbital elements for the selected objects
- `failed_horizons_ids.txt` — optional log of designations that failed Horizons queries

After running `split_rebound_cli.py`, the outputs are a sequence of directories such as:

- `split_1/`
- `split_2/`
- `split_3/`

Each split directory contains:

- a chunk of `ast_VRI.txt`
- the matching chunk of `Particles.el`
- a copied simulation template script, such as `Internal_MMR_ML_FAIR.py`

## Typical workflow

### Step 1: generate initial-condition files

```bash
python select_ast_earth_cli.py \
  --input-csv sbdb_query_results.csv \
  --a-max 0.98 \
  --epoch 2460800.5 \
  --location 500@10
```

This command:

- reads the SBDB query table
- keeps objects with `a < 0.98`
- sorts by `a`
- writes a filtered asteroid list to `ast_VRI.txt`
- queries JPL Horizons for each selected object
- writes orbital elements to `Particles.el`

### Step 2: split for simulation batches

```bash
python split_rebound_cli.py \
  --ast-file ast_VRI.txt \
  --particles-file Particles.el \
  --template-script Internal_MMR_ML_FAIR.py \
  --lines-per-batch 100
```

This command creates numbered batch directories and copies the simulation template into each one.

## Input data expectations

### SBDB CSV table

`select_ast_earth_cli.py` expects an input CSV with at least these columns:

- `a` — semi-major axis
- `pdes` — designation used to query JPL Horizons

Additional columns can also be written to the asteroid list output if requested with `--ast-columns`.

### Simulation template script

`split_rebound_cli.py` requires a Python template script that will be copied into every split directory. In the provided wrapper script this is:

- `Internal_MMR_ML_FAIR.py`

## Installation

Create and activate a virtual environment, then install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Notes for repository organization

For a publication-ready repository branch, a simple structure is:

```text
.
├── README.md
├── INPUTS_AND_FLAGS.md
├── requirements.txt
├── select_ast_earth_cli.py
├── split_rebound_cli.py
├── run_select_earth
├── run_split_rebound_cli
├── sbdb_query_results.csv          # optional example input, often excluded from Git
└── Internal_MMR_ML_FAIR.py         # simulation template copied into split directories
```

## Reproducibility notes

- Horizons queries depend on external service availability.
- The exact contents of `Particles.el` depend on the chosen epoch and location.
- Large input tables may take time because the current implementation performs one Horizons query per object.
- For reproducible batch runs, keep the same CSV input, filtering threshold, epoch, and template script.

## Citation / provenance

This README is based on the attached CLI scripts and wrapper commands for the branch setup. fileciteturn1file1 fileciteturn1file0
