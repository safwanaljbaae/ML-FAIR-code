# External MMR ML-FAIR

Machine-learning-assisted FAIR workflow for identifying **external mean-motion resonances (MMRs)** from asteroid integrations.

This branch documents the current behavior of the provided script `External_MMR_ML_FAIR.py` and its companion input files.

## What the code does

For each asteroid in the input sample, the script:

1. reads asteroid identifiers and reference semimajor axes from `ast_VRI.txt`;
2. reads asteroid orbital elements from `Particles.el`;
3. transforms orbital elements from the ecliptic of J2000 to the invariable plane;
4. queries planetary orbital elements from JPL Horizons;
5. integrates the Solar System + asteroid with `REBOUND`;
6. computes the external-resonance FAIR angle
   \[
   \sigma = (\Omega + \omega + M) - (\Omega_c + \omega_c + M_c)
   \]
   in the invariable plane, i.e. the asteroid mean longitude minus the selected planet mean longitude;
7. analyzes the angle distributions using KDE peak counting and OPTICS clustering;
8. flags ambiguous cases when the KDE and OPTICS counts disagree;
9. interprets the accepted FAIR counts for the external-resonance mapping
   - `p = nsigma`
   - `q = nM`
   - `p+q = nsigma + nM`
10. computes the theoretical external resonance semimajor axis
    \[
    a_{\mathrm{res}} = a_{\mathrm{planet}} \left(\frac{p+q}{p}\right)^{2/3}
    \]
11. writes a summary CSV for resonant candidates with `Delta_a < 0.4 AU`.

## Repository contents

Suggested minimum layout for this branch:

```text
external-mmr/
├── External_MMR_ML_FAIR.py
├── ast_VRI.txt
├── Particles.el
├── requirements.txt
├── INPUTS_AND_FLAGS.md
├── README.md
└── examples/
    ├── resonant_summary_mercury.csv
    ├── resonant_summary_venus.csv
    └── ...
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
0.8445 "2022 NX"
0.8452 467336
0.8508 480820
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
0.8445000000 0.1234567890 5.4321 120.1234 210.5678 45.6789
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
python External_MMR_ML_FAIR.py --planet mercury
python External_MMR_ML_FAIR.py --planet venus
python External_MMR_ML_FAIR.py --planet earth
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

If no `--planet` flag is provided, the script prompts:

```text
Choose a planet (mercury, venus, earth, mars, jupiter, saturn, uranus, neptune):
```

## Outputs

For each asteroid, the script produces diagnostic plots such as:

- `res_forgacs_<planet>_<asteroid>.png`
- `res_forgacs_<planet>_<asteroid>_FLAGGED.png`
- `res_forgacs_<planet>_<asteroid>_FAIR_FLAGGED.png`

For successful resonant identifications, it also writes:

```text
resonant_summary_<planet>.csv
```

This CSV includes, at minimum:
- asteroid identifier;
- selected planet;
- inferred resonance integers `p` and `p+q` (`pq` in the script);
- input semimajor axis;
- theoretical resonance semimajor axis for the selected planet;
- `Delta_a`;
- whether the case was flagged.

The summary file is only written for asteroids that satisfy both:
- finite resonance integers with `p > 0` and `q >= 0`; and
- `Delta_a < 0.4 AU`.

## Manual intervention for flagged asteroids

Some asteroids are marked as **FLAGGED** when the two counting methods disagree:

- KDE peak count versus OPTICS cluster count in `sigma @ M ≈ 0`, and/or
- KDE peak count versus OPTICS cluster count in `M @ sigma ≈ 0`.

When that happens, the script:

1. saves diagnostic figures;
2. opens the FAIR-only figure if possible;
3. pauses and asks the user to enter:
   - `nsigma`
   - `nM`

Prompt shown by the code:

```text
Enter manual nsigma (integer):
Enter manual nM (integer):
```

For the external code, those manual values are used as:

```text
p = nsigma
q = nM
pq = p + q
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
- multiprocessing is currently hard-coded to `max_workers=6`;
- output filenames do not include the word `external`, so they can collide with other runs if you reuse the same directory structure;
- the code uses fixed filtering and clustering thresholds in the script body rather than command-line options.

## Recommended next cleanup steps

1. move input and output paths into command-line arguments;
2. add a non-interactive mode for flagged asteroids;
3. add a small test dataset;
4. split the monolithic script into modules (`io`, `dynamics`, `fair`, `plotting`, `cli`);
5. add a license and citation information.

## Acknowledgment

This README reflects the current behavior of the provided research script and is intended as a starting point for the `external-mmr` branch.
