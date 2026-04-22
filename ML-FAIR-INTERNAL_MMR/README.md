# Inner MMR ML-FAIR

Machine-learning-assisted FAIR workflow for identifying **internal mean-motion resonances (MMRs)** from asteroid integrations.

This branch contains the first working version of the inner-resonance pipeline based on the provided script `Internal_MMR_ML_FAIR.py` and its companion input files.

## What the code does

For each asteroid in the input sample, the script:

1. reads asteroid identifiers and reference semimajor axes from `ast_VRI.txt`;
2. reads asteroid orbital elements from `Particles.el`;
3. transforms orbital elements from the ecliptic of J2000 to the invariable plane;
4. queries planetary orbital elements from JPL Horizons;
5. integrates the Solar System + asteroid with `REBOUND`;
6. computes the inner-resonance FAIR angle
   \[
   \sigma = (\Omega_c + \omega_c + M_c) - (\Omega + \omega + M)
   \]
   in the invariable plane;
7. analyzes the angle distributions using KDE peak counting and OPTICS clustering;
8. flags ambiguous cases when the KDE and OPTICS counts disagree;
9. converts the accepted FAIR counts into the internal resonance mapping
   - `nsigma = p + q`
   - `nM = q`
   - `p = nsigma - nM`
10. writes a summary CSV for resonant candidates with `Delta_a < 0.4 AU`.

## Repository contents

Suggested minimum layout for this branch:

```text
inner-mmr/
├── Internal_MMR_ML_FAIR.py
├── ast_VRI.txt
├── Particles.el
├── requirements.txt
├── INPUTS_AND_FLAGS.md
├── README.md
└── examples/
    ├── resonant_summary_internal_earth.csv
    └── resonant_summary_internal_mars.csv
```

## Dependencies

Install the Python dependencies with:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Required input files

The script expects the following files in the **working directory** unless the code is edited.

### 1. `ast_VRI.txt`
Plain-text asteroid list. Each non-comment line must contain:

```text
<a_input_AU> <asteroid_identifier>
```

Examples:

```text
0.8259 337248
0.8261 "2020 NC"
0.8269 "2004 FU162"
```

Notes:
- the first column is the reference/input semimajor axis in AU;
- the second field is the asteroid name or designation;
- quoted designations are accepted;
- spaces are stripped from the stored asteroid identifier inside the script.

### 2. `Particles.el`
Plain-text table of asteroid orbital elements, one asteroid per line, in this order:

```text
a e i Omega omega M
```

where angles are expected in degrees and `a` is expected in AU.

Example:

```text
0.825890258026559 0.5512851706493299 19.64965686446971 177.859822300090002 354.412106477580778 118.906332033273600
```

Important:
- the order of objects in `Particles.el` must match the order of asteroid IDs in `ast_VRI.txt`;
- the script pairs rows by index, not by name.

### 3. Optional `real`
Optional ground-truth file used only for evaluation metrics.
Each line is interpreted as:

```text
<asteroid_name> <n_clusters_sigma> <n_clusters_M>
```

If the file is absent, the code still runs, but accuracy metrics are skipped.

## Running the code

You can select the planet interactively or pass it on the command line.

Examples:

```bash
python Internal_MMR_ML_FAIR.py --planet earth
python Internal_MMR_ML_FAIR.py --planet mars
```

Valid planets:

- mercury
- venus
- earth
- mars
- jupiter
- saturn
- uranus
- neptune

## Outputs

For each asteroid, the script produces diagnostic plots such as:

- `res_forgacs_internal_<planet>_<asteroid>.png`
- `res_forgacs_internal_<planet>_<asteroid>_FLAGGED.png`
- `res_forgacs_internal_<planet>_<asteroid>_FAIR_FLAGGED.png`

For successful resonant identifications, it also writes:

```text
resonant_summary_internal_<planet>.csv
```

This CSV includes, at minimum:
- asteroid identifier;
- selected planet;
- inferred resonance integers `p` and `p+q` (`pq` in the script);
- input semimajor axis;
- theoretical resonance semimajor axis;
- `Delta_a`;
- whether the case was flagged.

## Manual intervention for flagged asteroids

Some asteroids are marked as **FLAGGED** when the two counting methods disagree:

- KDE peak count versus OPTICS cluster count in `sigma @ M ≈ 0`, and/or
- KDE peak count versus OPTICS cluster count in `M @ sigma ≈ 0`.

When that happens, the script:

1. saves diagnostic figures;
2. opens the flagged figure if possible;
3. pauses and asks the user to enter:
   - `nsigma` = the adopted number of vertical FAIR intersections (`p + q`);
   - `nM` = the adopted number of horizontal FAIR intersections (`q`).

Prompt shown by the code:

```text
Enter manual nsigma (integer):
Enter manual nM (integer):
```

This means the current workflow is **semi-automatic**, not fully unattended.
For batch use on a cluster or CI system, you may want to replace this prompt with:

- a review table prepared in advance;
- a YAML/CSV override file for flagged asteroids; or
- a rule that skips flagged objects until manual inspection is completed.

## Current limitations

This branch is already enough to start, but a few issues should be kept in mind:

- the plotting backend is hard-coded as `TkAgg`, which may fail on headless servers;
- the script assumes fixed filenames in the working directory;
- asteroid matching between `ast_VRI.txt` and `Particles.el` is positional only;
- the JPL Horizons query requires internet access at runtime;
- multiprocessing is currently hard-coded to `max_workers=6`.

## Recommended next cleanup steps

1. move input and output paths into command-line arguments;
2. add a non-interactive mode for flagged asteroids;
3. add a small test dataset;
4. split the monolithic script into modules (`io`, `dynamics`, `fair`, `plotting`, `cli`);
5. add a license and citation information.

## Acknowledgment

This README reflects the current behavior of the provided research script and is intended as a starting point for the `inner-mmr` branch.
