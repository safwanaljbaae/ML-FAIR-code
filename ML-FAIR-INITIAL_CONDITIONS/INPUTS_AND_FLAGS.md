# Inputs and Flags

This document summarizes the required inputs, generated outputs, and command-line flags for the initial-condition preparation scripts in this branch.

## 1. `select_ast_earth_cli.py`

### Purpose

Filter an SBDB-style asteroid catalog and build the two main input files used by the simulation workflow:

- `ast_VRI.txt`
- `Particles.el`

### Required external input

- **SBDB CSV file**
  - provided with `--input-csv`
  - must contain at least:
    - `a`
    - `pdes`

### Main outputs

- **Asteroid list file**
  - default: `ast_VRI.txt`
  - controlled by `--ast-output`
- **Orbital elements file**
  - default: `Particles.el`
  - controlled by `--particles-output`
- **Failed Horizons ID log**
  - default: `failed_horizons_ids.txt`
  - controlled by `--failed-ids-output`

### Flags

#### `--input-csv`
Path to the SBDB-style CSV input table.

Default:
```text
sbdb_query_results.csv
```

#### `--a-max`
Upper threshold for the semi-major axis filter. The script keeps rows satisfying:

```text
a < a_max
```

Default:
```text
0.98
```

#### `--sort-by`
Column used to sort the filtered table before writing outputs.

Default:
```text
a
```

#### `--designation-column`
Column used as the Horizons object identifier.

Default:
```text
pdes
```

#### `--ast-columns`
Comma-separated list of columns written to the asteroid list file.

Default:
```text
a,pdes
```

#### `--ast-output`
Output path for the filtered asteroid list.

Default:
```text
ast_VRI.txt
```

#### `--particles-output`
Output path for the orbital-elements file.

Default:
```text
Particles.el
```

#### `--epoch`
Julian Date epoch used for JPL Horizons element queries.

Default:
```text
2460800.5
```

#### `--location`
Location / center code passed to JPL Horizons.

Default:
```text
500@10
```

#### `--limit`
Optional maximum number of filtered asteroids to process.

Default:
```text
None
```

This is useful for:

- quick tests
- debugging
- small demonstration runs

#### `--sleep`
Seconds to wait between Horizons queries.

Default:
```text
0.0
```

Useful when you want to reduce API pressure during large runs.

#### `--failed-ids-output`
Output file where failed Horizons identifiers are written.

Default:
```text
failed_horizons_ids.txt
```

### Example

```bash
python select_ast_earth_cli.py \
  --input-csv sbdb_query_results.csv \
  --a-max 0.98 \
  --sort-by a \
  --designation-column pdes \
  --ast-columns a,pdes \
  --ast-output ast_VRI.txt \
  --particles-output Particles.el \
  --epoch 2460800.5 \
  --location 500@10
```

## 2. `split_rebound_cli.py`

### Purpose

Split `ast_VRI.txt` and `Particles.el` into a set of numbered directories for batched simulation runs, and copy a template simulation script into each directory.

### Required inputs

- asteroid list file (`--ast-file`)
- orbital elements file (`--particles-file`)
- simulation template script (`--template-script`)

### Main outputs

- numbered directories such as `split_1`, `split_2`, `split_3`, ...
- each directory contains:
  - one chunk of the asteroid list
  - one chunk of the particles file
  - one copied simulation script

### Flags

#### `--ast-file`
Path to the asteroid list file.

Default:
```text
ast_VRI.txt
```

#### `--particles-file`
Path to the orbital elements file.

Default:
```text
Particles.el
```

#### `--template-script`
Python script copied into every generated batch directory.

Default:
```text
Earth_inner_ensemble.py
```

In the provided branch wrapper, this is set to:

```text
Internal_MMR_ML_FAIR.py
```

#### `--lines-per-batch`
Maximum number of asteroid / particle rows in each output directory.

Default:
```text
100
```

#### `--output-prefix`
Prefix used to name generated directories.

Default:
```text
split_
```

Example output directory names:

```text
split_1
split_2
split_3
```

#### `--copied-script-name`
Filename to use for the copied template script inside each split directory.

Default:
```text
Earth_inner_ensemble.py
```

#### `--strict`
If supplied, the script fails when `--ast-file` and `--particles-file` do not contain the same number of lines.

Default behavior:
```text
Disabled
```

#### `--overwrite`
If supplied, existing output directories are deleted and recreated.

Default behavior:
```text
Disabled
```

### Example

```bash
python split_rebound_cli.py \
  --ast-file ast_VRI.txt \
  --particles-file Particles.el \
  --template-script Internal_MMR_ML_FAIR.py \
  --lines-per-batch 100 \
  --overwrite
```

## 3. Wrapper scripts included in the branch

### `run_select_earth`
A short convenience script that calls:

```bash
python select_ast_earth_cli.py \
  --input-csv sbdb_query_results.csv \
  --a-max 0.98 \
  --epoch 2460800.5 \
  --location 500@10
```

### `run_split_rebound_cli`
A short convenience script that calls:

```bash
python split_rebound_cli.py \
  --ast-file ast_VRI.txt \
  --particles-file Particles.el \
  --template-script Internal_MMR_ML_FAIR.py \
  --lines-per-batch 100
```

## 4. Practical recommendations for the GitHub branch

- Commit the Python source files and documentation.
- Usually do **not** commit generated outputs such as:
  - `ast_VRI.txt`
  - `Particles.el`
  - `failed_horizons_ids.txt`
  - `split_*/`
- Commit example wrapper scripts if you want users to reproduce your standard workflow quickly.
- Consider excluding large raw input tables unless you specifically want to distribute an example dataset.

## Source basis

This flag summary is derived from the attached CLI scripts and wrapper commands. fileciteturn1file1 fileciteturn1file0
