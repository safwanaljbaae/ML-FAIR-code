# Notes on input files and flagged asteroids

## Required input files

The current version of `Internal_MMR_ML_FAIR.py` expects local files with fixed names in the working directory.

### `ast_VRI.txt`
This file defines the asteroid sample and the reference semimajor axis used later in the resonance summary.

Expected format per non-comment line:

```text
<a_input_AU> <asteroid_identifier>
```

Examples:

```text
0.8259 337248
0.8261 "2020 NC"
0.8269 "2004 FU162"
```

Behavior in the code:
- blank lines and lines beginning with `#` are ignored;
- the first token is parsed as a float and stored as `a_input`;
- the remainder of the line is taken as the asteroid identifier;
- quotes are stripped;
- embedded spaces are removed from the stored identifier.

### `Particles.el`
This file supplies the asteroid orbital elements used by the N-body simulation.

Expected format per line:

```text
a e i Omega omega M
```

Units expected by the code:
- `a` in AU
- `e` dimensionless
- angular elements in degrees

Example:

```text
0.825890258026559 0.5512851706493299 19.64965686446971 177.859822300090002 354.412106477580778 118.906332033273600
```

Important coupling rule:
- `Particles.el` and `ast_VRI.txt` must describe the same asteroid list in the same row order.
- The script does **not** cross-match by name.

### Optional `real`
Optional evaluation file for comparing automated counts against trusted labels.

Expected format:

```text
<asteroid_name> <n_clusters_sigma> <n_clusters_M>
```

If omitted, the analysis still runs and only the accuracy metrics section is skipped.

## What is a flagged asteroid?

An asteroid is flagged when the automated counting methods do not agree.

The script compares:
- the number of KDE peaks in `sigma @ M ≈ 0` with the number of OPTICS clusters in the same distribution;
- the number of KDE peaks in `M @ sigma ≈ 0` with the number of OPTICS clusters in that distribution.

If either comparison disagrees, the asteroid is marked `FLAGGED`.

## What the user must do for flagged asteroids

For each flagged asteroid, the script prints diagnostic information, writes plot files, and then asks the user to inspect the case manually.

The user must supply two integers:

- `nsigma`: interpreted in the script as `p + q`
- `nM`: interpreted in the script as `q`

The final inner resonance order is then computed as:

```text
p = nsigma - nM
pq = nsigma
```

## Practical recommendation for repository users

Please make it explicit in the branch documentation that the current version is **semi-automatic**.
A complete run may stop and wait for human input.

For reproducible production runs, consider one of these improvements:

1. maintain a manual override file for flagged asteroids;
2. store reviewer decisions in CSV form;
3. add a `--non-interactive` mode that skips unresolved flagged cases;
4. log all manual choices for provenance.
