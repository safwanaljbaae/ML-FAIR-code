import argparse
import concurrent.futures
import os

import astropy.table
import astropy.units as u
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rebound
from angles import normalize
from astroquery.jplhorizons import Horizons
from scipy.ndimage import maximum_filter1d
from scipy.signal import peak_prominences
from scipy.stats import gaussian_kde
from sklearn.cluster import OPTICS
from sklearn.metrics import accuracy_score

# -----------------------------------------------------------
# User options
# -----------------------------------------------------------
# You can choose the planet interactively, or pass it via:
#   python Internal_MMR_ML_FAIR.py --planet mercury
# Valid options: mercury, venus, earth, mars, jupiter, saturn, uranus, neptune

# -----------------------------------------------------------
# rotation_matrix import (compatible across astropy versions)
# -----------------------------------------------------------
try:
    from astropy.coordinates.matrix_utilities import rotation_matrix
except ImportError:
    from astropy.coordinates.angles import rotation_matrix

# ============================================================
# Invariable Plane Transformation Utilities
# ============================================================

I_INV = 1.5833333 * u.deg
OMEGA_INV = 107.0 * u.deg


def to_invariable_plane(a, e, inc, Omega, omega):
    """Transform orbital elements from ecliptic of J2000 to invariable plane."""
    i_rad = np.radians(inc)
    Omega_rad = np.radians(Omega)

    Rz = rotation_matrix(-OMEGA_INV, "z")
    Rx = rotation_matrix(-I_INV, "x")
    R = Rx @ Rz

    h_vec = np.array([
        np.sin(i_rad) * np.sin(Omega_rad),
        -np.sin(i_rad) * np.cos(Omega_rad),
        np.cos(i_rad),
    ])
    h_vec_new = R @ h_vec

    i_new = np.degrees(np.arccos(h_vec_new[2]))
    Omega_new = np.degrees(np.arctan2(h_vec_new[0], -h_vec_new[1])) % 360.0
    omega_new = (omega + Omega - Omega_new) % 360.0
    return a, e, i_new, Omega_new, omega_new


# ============================================================
# Read asteroid names + input semi-major axis
# Format per line: 0.8392 "2010 MB"
# ============================================================

input_file = "ast_VRI.txt"
ast_id = []
ast_a_input = []

with open(input_file, "r") as fin:
    for raw in fin:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split(maxsplit=1)
        if len(parts) < 2:
            continue

        try:
            a_input = float(parts[0])
        except ValueError:
            continue

        name_clean = parts[1].strip().strip('"').replace(" ", "")
        ast_a_input.append(a_input)
        ast_id.append(name_clean)

print(f"Loaded {len(ast_id)} asteroids from {input_file}")


# ============================================================
# Read real (ground-truth) peak/cluster counts for evaluation
# ============================================================

def load_real_results(path="real"):
    """Load ground-truth cluster counts from a file.

    Parameters
    ----------
    path : str, optional
        Path to the file containing ground-truth data. Default is "real".

    Returns
    -------
    dict
        Dictionary mapping asteroid names to their ground-truth cluster counts.
        Each entry contains keys 'n_clusters_sigma' and 'n_clusters_M'.
    """
    real_map = {}
    if not os.path.exists(path):
        return real_map
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            ast_name = parts[0].strip().replace(" ", "")
            try:
                n_sig = int(float(parts[1]))
                n_m = int(float(parts[2]))
            except ValueError:
                continue
            real_map[ast_name] = {"n_clusters_sigma": n_sig, "n_clusters_M": n_m}
    return real_map


# ============================================================
# Read asteroid orbital elements and transform to invariable plane
# ============================================================

def ler_elementos_particles(arquivo="Particles.el"):
    """Read asteroid orbital elements from a file and transform to invariable plane.

    Parameters
    ----------
    arquivo : str, optional
        Path to the file containing orbital elements. Default is "Particles.el".

    Returns
    -------
    list of tuple
        List of tuples containing (a, e, inc_rad, Omega_rad, omega_rad, M_rad)
        for each asteroid, with angles in radians and elements transformed to the
        invariable plane.
    """
    with open(arquivo, "r") as f:
        linhas = f.readlines()
    elementos = []
    for linha in linhas:
        partes = linha.strip().split()
        a, e, inc, Omega, omega, M = map(float, partes[:6])
        a, e, inc_new, Omega_new, omega_new = to_invariable_plane(a, e, inc, Omega, omega)
        elementos.append((
            a,
            e,
            np.radians(inc_new),
            np.radians(Omega_new),
            np.radians(omega_new),
            np.radians(M),
        ))
    return elementos


# ============================================================
# Planets data
# ============================================================

bodies = ["199", "299", "399", "301", "499", "599", "699", "799", "899"]
all_results = astropy.table.vstack([
    Horizons(id=i, location="500@10", epochs=2460800.5).elements()
    for i in bodies
])

for row in all_results:
    a, e = row["a"], row["e"]
    inc, Omega, omega = row["incl"], row["Omega"], row["w"]
    a, e, inc_new, Omega_new, omega_new = to_invariable_plane(a, e, inc, Omega, omega)
    row["incl"] = inc_new
    row["Omega"] = Omega_new
    row["w"] = omega_new

# ============================================================
# Planetary masses (solar masses)
# ============================================================

massas = [
    1.6601208254808336e-07,  # Mercury
    2.447838287784771e-06,   # Venus
    3.0034896154502038e-06,  # Earth
    3.694303350091508e-08,   # Moon
    3.2271559174983593e-07,  # Mars
    9.545942479871165e-04,   # Jupiter
    2.8581500081698117e-04,  # Saturn
    4.365793681773934e-05,   # Uranus
    5.1503084159155005e-05,  # Neptune
]

# ============================================================
# Planet configuration
# sim_index = particle index in REBOUND after adding the asteroid
# all_results_index = row index in all_results
# ============================================================

PLANET_INFO = {
    "mercury": {"sim_index": 1, "all_results_index": 0, "label": "Mercury"},
    "venus":   {"sim_index": 2, "all_results_index": 1, "label": "Venus"},
    "earth":   {"sim_index": 3, "all_results_index": 2, "label": "Earth"},
    "mars":    {"sim_index": 5, "all_results_index": 4, "label": "Mars"},
    "jupiter": {"sim_index": 7, "all_results_index": 5, "label": "Jupiter"},
    "saturn":  {"sim_index": 8, "all_results_index": 6, "label": "Saturn"},
    "uranus":  {"sim_index": 9, "all_results_index": 7, "label": "Uranus"},
    "neptune": {"sim_index": 10, "all_results_index": 8, "label": "Neptune"},
}


# ============================================================
# Simulation time settings
# ============================================================

total_days = int(4.1e5)
total_years = total_days / 365.25
dt = 0.1095
times = np.arange(0.0, total_years + dt, dt)
N = len(times)


# ============================================================
# Resonant angle in invariable plane (internal case)
# ============================================================

def compute_sigma_elements_invariable(planet, asteroid, primary):
    """Compute the resonant angle sigma in the invariable plane for internal resonances.

    Sigma is defined as the difference between the planet's longitude and the asteroid's
    longitude: sigma = (Omega + omega + M)_planet - (Omega + omega + M)_asteroid.
    For internal resonances, the asteroid is inside the planet's orbit.

    Parameters
    ----------
    planet : rebound.Particle
        The planet particle.
    asteroid : rebound.Particle
        The asteroid particle.
    primary : rebound.Particle
        The primary body (Sun).

    Returns
    -------
    float
        The resonant angle sigma in radians, normalized to [0, 2*pi).
    """
    orb_p = planet.orbit(primary=primary)
    orb_a = asteroid.orbit(primary=primary)

    lam_p = (orb_p.Omega + orb_p.omega + orb_p.M) % (2 * np.pi)
    lam_a = (orb_a.Omega + orb_a.omega + orb_a.M) % (2 * np.pi)

    return (lam_p - lam_a) % (2 * np.pi)


# ============================================================
# Circular union coverage
# ============================================================

def circular_union_coverage(angles_deg, halfwidth_deg, period=360.0):
    """Calculate the fractional coverage of a circular domain by intervals around given angles.

    For each angle in the input list, creates an interval of halfwidth_deg on either side,
    then computes the union of all intervals and returns the fraction of the circle covered.

    Parameters
    ----------
    angles_deg : array-like
        List of center angles in degrees.
    halfwidth_deg : float
        Half-width of each interval in degrees.
    period : float, optional
        Total period of the circle in degrees. Default is 360.0.

    Returns
    -------
    tuple
        (fraction_covered, merged_intervals, gaps)
        - fraction_covered: Fraction of the circle covered by the union of intervals.
        - merged_intervals: List of (start, end) tuples for merged intervals.
        - gaps: List of (start, end) tuples for gaps in coverage.
    """
    ang = np.asarray(angles_deg, dtype=float)
    if ang.size == 0:
        return 0.0, [], [(0.0, period)]

    w = float(halfwidth_deg)
    ang = np.mod(ang, period)
    intervals = []

    for a in ang:
        lo = a - w
        hi = a + w
        if lo < 0:
            intervals.append((lo + period, period))
            intervals.append((0.0, hi))
        elif hi >= period:
            intervals.append((lo, period))
            intervals.append((0.0, hi - period))
        else:
            intervals.append((lo, hi))

    intervals.sort(key=lambda x: x[0])
    merged = []
    cur_lo, cur_hi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            merged.append((cur_lo, cur_hi))
            cur_lo, cur_hi = lo, hi
    merged.append((cur_lo, cur_hi))

    covered_len = sum(hi - lo for lo, hi in merged)
    frac = covered_len / period

    gaps = []
    prev = 0.0
    for lo, hi in merged:
        if lo > prev:
            gaps.append((prev, lo))
        prev = hi
    if prev < period:
        gaps.append((prev, period))

    return frac, merged, gaps


# ============================================================
# OPTICS for circular angles
# ============================================================

def optics_circular_angles(angles_deg, max_eps_deg=10.0, min_samples=5):
    """Cluster circular angles using the OPTICS algorithm.

    Converts angles to 2D Cartesian coordinates on the unit circle and applies
    OPTICS clustering to identify groups of nearby angles, accounting for
    circular periodicity.

    Parameters
    ----------
    angles_deg : array-like
        List of angles in degrees.
    max_eps_deg : float, optional
        Maximum epsilon distance for clustering in degrees. Default is 10.0.
    min_samples : int, optional
        Minimum number of points to form a cluster. Default is 5.

    Returns
    -------
    tuple
        (labels, n_clusters)
        - labels: Array of cluster labels for each angle.
        - n_clusters: Number of clusters found.
    """
    angles_deg = np.asarray(angles_deg, dtype=float)
    if angles_deg.size == 0:
        return np.array([], dtype=int), 0

    angles_rad = np.deg2rad(angles_deg)
    X = np.column_stack((np.cos(angles_rad), np.sin(angles_rad)))
    max_eps = 2.0 * np.sin(np.deg2rad(max_eps_deg) / 2.0)

    optics = OPTICS(
        min_samples=min_samples,
        max_eps=max_eps,
        cluster_method="dbscan",
        eps=max_eps,
    )
    labels = optics.fit_predict(X)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    return labels, n_clusters


# ============================================================
# Circular KDE + peak counting
# ============================================================

def _cluster_boolean_runs_circular(is_max):
    """Identify representative indices for boolean runs in a circular array.

    Helper function for count_peaks_circular. Finds the center point of each
    contiguous run of True values in a circular boolean array.

    Parameters
    ----------
    is_max : array-like
        Boolean array indicating local maxima.

    Returns
    -------
    np.ndarray
        Array of representative indices for each run of True values.
    """
    n = len(is_max)
    idx = np.flatnonzero(is_max)
    if idx.size == 0:
        return np.array([], dtype=int)

    splits = np.where(np.diff(idx) > 1)[0] + 1
    groups = np.split(idx, splits)

    if groups and groups[0][0] == 0 and groups[-1][-1] == n - 1:
        groups[0] = np.r_[groups[-1], groups[0]]
        groups = groups[:-1]

    reps = np.array([g[len(g) // 2] % n for g in groups], dtype=int)
    return reps


def count_peaks_circular(data_deg, bandwidth_deg=10.0, prom_frac_range=0.2, n_grid=2000):
    """Count peaks in a circular distribution using kernel density estimation.

    Applies KDE to the angular data, then finds peaks by identifying local maxima
    with prominence above a threshold. Handles circular periodicity by extending
    the data with wrapped copies.

    Parameters
    ----------
    data_deg : array-like
        List of angles in degrees.
    bandwidth_deg : float, optional
        Bandwidth for the KDE in degrees. Default is 10.0.
    prom_frac_range : float, optional
        Minimum prominence as fraction of KDE dynamic range. Default is 0.2.
    n_grid : int, optional
        Number of grid points for KDE evaluation. Default is 2000.

    Returns
    -------
    tuple
        (n_peaks, peak_positions, grid_deg, kde_values)
        - n_peaks: Number of peaks detected.
        - peak_positions: Array of peak positions in degrees.
        - grid_deg: Grid of angle values in degrees.
        - kde_values: KDE values at each grid point.
    """
    data_deg = np.asarray(data_deg, dtype=float)
    if data_deg.size < 5:
        return 0, np.array([]), np.array([]), np.array([])

    data_ext = np.concatenate([data_deg, data_deg + 360.0, data_deg - 360.0])
    kde = gaussian_kde(np.deg2rad(data_ext), bw_method=np.deg2rad(bandwidth_deg))

    grid_deg = np.linspace(0.0, 360.0, n_grid, endpoint=False)
    kde_vals = kde(np.deg2rad(grid_deg))

    dyn = float(kde_vals.max() - kde_vals.min())
    if dyn <= 0:
        return 0, np.array([]), grid_deg, kde_vals
    prom_abs = prom_frac_range * dyn

    step = 360.0 / n_grid
    radius = max(1, int(np.ceil(bandwidth_deg / step)))
    win = 2 * radius + 1

    y_max = maximum_filter1d(kde_vals, size=win, mode="wrap")
    is_max = kde_vals >= y_max - 1e-12
    cand = _cluster_boolean_runs_circular(is_max)

    y_ext = np.r_[kde_vals, kde_vals, kde_vals]
    cand_ext = cand + n_grid
    prom = peak_prominences(y_ext, cand_ext)[0]

    keep = prom >= prom_abs
    peaks = cand[keep]

    peak_pos = (peaks * step) % 360.0
    peak_pos.sort()
    return len(peaks), peak_pos, grid_deg, kde_vals


# ============================================================
# Image display helper
# ============================================================

def show_image_matplotlib(path: str):
    """Display an image file using matplotlib.

    Parameters
    ----------
    path : str
        Path to the image file to display.

    Returns
    -------
    matplotlib.figure.Figure or None
        The created figure, or None if display failed.
    """
    try:
        img = plt.imread(path)
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(img)
        ax.axis("off")
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.1)
        return fig
    except Exception as e:
        print(f"Could not display image {path}: {e}")
        return None


# ============================================================
# Inner resonance location using Kepler's third law
# ============================================================

def resonance_a_inner(a_planet_au: float, p: int, pq: int) -> float:
    """Calculate the semi-major axis of an inner p:q mean motion resonance.

    Uses Kepler's third law to compute the resonance location: a_res = a_planet * (p/pq)^(2/3).
    For inner resonances, the asteroid is inside the planet's orbit.

    Parameters
    ----------
    a_planet_au : float
        Semi-major axis of the planet in AU.
    p : int
        Resonance integer p (non-zero).
    pq : int
        Resonance integer pq (non-zero).

    Returns
    -------
    float
        Semi-major axis of the resonance in AU, or NaN if inputs are invalid.
    """
    if p <= 0 or pq <= 0:
        return float("nan")
    return a_planet_au * (p / pq) ** (2.0 / 3.0)


# ============================================================
# Planet selection
# ============================================================

def get_planet_name():
    """Get the planet name from command-line arguments or interactive input.

    Returns
    -------
    str
        Lowercase planet name (e.g., 'mercury', 'venus', 'earth', etc.).

    Raises
    ------
    ValueError
        If an invalid planet name is provided.
    """
    parser = argparse.ArgumentParser(
        description="Compute internal MMR locations with ML-FAIR for a selected planet."
    )
    parser.add_argument(
        "--planet",
        type=str,
        help="Planet name: mercury, venus, earth, mars, jupiter, saturn, uranus, neptune",
    )
    args = parser.parse_args()

    if args.planet:
        planet_name = args.planet.strip().lower()
    else:
        valid = ", ".join(PLANET_INFO.keys())
        planet_name = input(f"Choose a planet ({valid}): ").strip().lower()

    if planet_name not in PLANET_INFO:
        raise ValueError(
            f"Invalid planet '{planet_name}'. Choose one of: {', '.join(PLANET_INFO.keys())}"
        )
    return planet_name


SELECTED_PLANET = get_planet_name()
PLANET_CFG = PLANET_INFO[SELECTED_PLANET]
print(f"Selected internal-resonance planet: {PLANET_CFG['label']}")


# ============================================================
# Single asteroid simulation
# ============================================================

def simular_asteroide(args):
    """Run N-body simulation for a single asteroid and analyze internal resonant behavior.

    Sets up a REBOUND simulation with the Sun, all planets, and the asteroid.
    Integrates the system and analyzes the resonant angle sigma and asteroid's
    mean anomaly M to detect internal mean motion resonances. Generates visualization
    plots and returns cluster counts.

    Parameters
    ----------
    args : tuple
        (idx, elementos_asteroide) where idx is the asteroid index and
        elementos_asteroide is a tuple of (a, e, inc, Omega, omega, M) in
        radians (already transformed to invariable plane).

    Returns
    -------
    tuple
        (idx, name, a_input, log_lines, filename, n_peaks_s, n_peaks_M,
         n_clusters_sigma, n_clusters_M, nsigma, nM, flagged, fair_filename)
    """
    idx, elementos_asteroide = args
    a_ast, e_ast, inc_ast, Omega_ast, omega_ast, M_ast = elementos_asteroide

    safe_idx = min(idx, len(ast_id) - 1)
    name = ast_id[safe_idx]
    a_input = ast_a_input[safe_idx]

    log_lines = [name]
    flagged = False
    nsigma = np.nan
    nM = np.nan

    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "bs"
    sim.add(m=0.9999999999950272)

    i = 0

    # Inner planets + Moon
    for row in all_results[0:5]:
        sim.add(
            m=massas[i],
            a=row["a"],
            e=row["e"],
            inc=np.radians(row["incl"]),
            Omega=np.radians(row["Omega"]),
            omega=np.radians(row["w"]),
            M=np.radians(row["M"]),
            primary=sim.particles[0],
        )
        i += 1

    # Asteroid
    sim.add(
        m=0,
        a=a_ast,
        e=e_ast,
        inc=inc_ast,
        Omega=Omega_ast,
        omega=omega_ast,
        M=M_ast,
        primary=sim.particles[0],
    )

    # Outer planets
    for row in all_results[5:]:
        sim.add(
            m=massas[i],
            a=row["a"],
            e=row["e"],
            inc=np.radians(row["incl"]),
            Omega=np.radians(row["Omega"]),
            omega=np.radians(row["w"]),
            M=np.radians(row["M"]),
            primary=sim.particles[0],
        )
        i += 1

    sim.move_to_com()

    time = np.zeros(N)
    sigma_time = np.zeros(N)
    M_Ast = np.zeros(N)

    for k, t in enumerate(times):
        sim.integrate(t)
        orb = sim.particles[6].orbit(primary=sim.particles[0])
        time[k] = sim.t
        M_Ast[k] = np.degrees(normalize(orb.M))
        sigma_rad = compute_sigma_elements_invariable(
            sim.particles[PLANET_CFG["sim_index"]],
            sim.particles[6],
            sim.particles[0],
        )
        sigma_time[k] = np.degrees(normalize(sigma_rad))

    M = np.mod(M_Ast, 360.0)
    sigma = np.mod(sigma_time, 360.0)

    tol = 5.0
    dens_cov = 0.945
    eps_deg = 15.0
    min_samples = 2
    bandwidth_deg = 2.5
    prom_frac_range = 0.40

    mask_M0 = (M <= tol) | (M >= 360.0 - tol)
    sigma_at_M0_list = sigma[mask_M0]

    mask_sigma0 = (sigma <= tol) | (sigma >= 360.0 - tol)
    M_at_sigma0_list = M[mask_sigma0]

    frac_sigma, _, _ = circular_union_coverage(sigma_at_M0_list, halfwidth_deg=2)
    frac_M, _, _ = circular_union_coverage(M_at_sigma0_list, halfwidth_deg=2)

    log_lines.append(f"σ@M≈0 coverage fraction: {frac_sigma}")
    log_lines.append(f"M@σ≈0 coverage fraction: {frac_M}")

    if (frac_sigma >= dens_cov) and (frac_M >= dens_cov):
        n_peaks_s = 0
        n_peaks_M = 0
        n_clusters_sigma = 0
        n_clusters_M = 0
        peak_pos_s = np.array([])
        peak_pos_M = np.array([])
        grid_s = grid_M = np.array([])
        kde_vals_s = kde_vals_M = np.array([])
        log_lines.append(
            f"High density coverage in both σ and M (>= {dens_cov}); skipping KDE and OPTICS."
        )
    else:
        n_peaks_s, peak_pos_s, grid_s, kde_vals_s = count_peaks_circular(
            sigma_at_M0_list,
            bandwidth_deg=bandwidth_deg,
            prom_frac_range=prom_frac_range,
        )
        log_lines.append(f"Number of peaks for σ at M≈0: {n_peaks_s}")

        n_peaks_M, peak_pos_M, grid_M, kde_vals_M = count_peaks_circular(
            M_at_sigma0_list,
            bandwidth_deg=5,
            prom_frac_range=0.40,
        )
        log_lines.append(f"Number of peaks for M at σ≈0: {n_peaks_M}")

        _, n_clusters_sigma = optics_circular_angles(
            sigma_at_M0_list,
            max_eps_deg=eps_deg,
            min_samples=min_samples,
        )
        log_lines.append(f"Number of clusters in σ at M≈0: {n_clusters_sigma}")

        _, n_clusters_M = optics_circular_angles(
            M_at_sigma0_list,
            max_eps_deg=eps_deg,
            min_samples=min_samples,
        )
        log_lines.append(f"Number of clusters in M at σ≈0: {n_clusters_M}")

        if (
            n_peaks_s == 1 and
            n_peaks_M == 1 and
            n_clusters_sigma == 1 and
            n_clusters_M == 1
        ):
            n_peaks_s = 0
            n_peaks_M = 0
            n_clusters_sigma = 0
            n_clusters_M = 0
            nsigma = 0
            nM = 0
            log_lines.append("Single peak/cluster in both σ and M → resetting all to 0.")

        if (
            (n_peaks_s > 0 and n_clusters_sigma > 0 and n_peaks_s != n_clusters_sigma) or
            (n_peaks_M > 0 and n_clusters_M > 0 and n_peaks_M != n_clusters_M)
        ):
            flagged = True
            log_lines.append("FLAGGED: KDE peaks and OPTICS clusters disagree (σ and/or M).")

    if (n_peaks_s == n_clusters_sigma) and (n_peaks_M == n_clusters_M):
        nsigma = int(n_peaks_s)
        nM = int(n_peaks_M)
        log_lines.append(f"Agreed counts → nsigma={nsigma}, nM={nM}")
    else:
        nsigma = np.nan
        nM = np.nan

    cvals = np.linspace(0.0, 1.0, len(time))
    fair_filename = None

    if flagged:
        fair_filename = f"res_forgacs_internal_{SELECTED_PLANET}_{name}_FAIR_FLAGGED.png"
        fig_fair, ax_fair = plt.subplots(1, 1, figsize=(5, 4))
        sc2 = ax_fair.scatter(M_Ast, sigma_time, c=cvals, cmap="coolwarm", marker=".", s=1)
        cbar2 = plt.colorbar(sc2, ax=ax_fair)
        cbar2.set_label("Normalized integration time", fontsize=12)
        ax_fair.set_xlabel("M_ast [deg]", fontsize=16)
        ax_fair.set_ylabel(r"$\sigma$ [deg]", fontsize=14)
        ax_fair.set_xlim([0, 360])
        ax_fair.set_ylim([0, 360])
        plt.tight_layout()
        plt.savefig(fair_filename)
        plt.close(fig_fair)

    fig, axs = plt.subplots(1, 5, figsize=(15, 4))

    sc = axs[0].scatter(M_Ast, sigma_time, c=cvals, cmap="coolwarm", marker=".", s=1)
    cbar = plt.colorbar(sc, ax=axs[0])
    cbar.set_label("Normalized integration time", fontsize=12)
    axs[0].set_xlabel("M_ast [deg]", fontsize=16)
    axs[0].set_ylabel(r"$\sigma = (\Omega_c+\omega_c+M_c)-(\Omega+\omega+M)$ [deg]", fontsize=14)
    axs[0].set_xlim([0, 360])
    axs[0].set_ylim([0, 360])

    axs[1].hist(sigma_at_M0_list, bins=30, color="tab:blue", alpha=0.7)
    axs[1].set_xlabel(r"$\sigma$ at $M=0$ [deg]")
    axs[1].set_ylabel("Count")
    axs[1].set_xlim([0, 360])

    axs[2].hist(M_at_sigma0_list, bins=30, color="tab:orange", alpha=0.7)
    axs[2].set_xlabel(r"$M$ at $\sigma=0$ [deg]")
    axs[2].set_ylabel("Count")
    axs[2].set_xlim([0, 360])

    axs[3].plot(grid_s, kde_vals_s, color="k")
    if peak_pos_s.size:
        axs[3].scatter(
            peak_pos_s,
            kde_vals_s[np.searchsorted(grid_s, peak_pos_s)],
            color="r",
            zorder=3,
        )
    axs[3].set_xlabel(r"$\sigma$ at $M \approx 0^\circ$ [deg]")
    axs[3].set_ylabel("KDE")

    axs[4].plot(grid_M, kde_vals_M, color="k")
    if peak_pos_M.size:
        axs[4].scatter(
            peak_pos_M,
            kde_vals_M[np.searchsorted(grid_M, peak_pos_M)],
            color="r",
            zorder=3,
        )
    axs[4].set_xlabel(r"$M$ at $\sigma \approx 0^\circ$ [deg]")
    axs[4].set_ylabel("KDE")

    plt.tight_layout()
    filename = f"res_forgacs_internal_{SELECTED_PLANET}_{name}{'_FLAGGED' if flagged else ''}.png"
    plt.savefig(filename)
    plt.close()

    return (
        idx,
        name,
        a_input,
        log_lines,
        filename,
        n_peaks_s,
        n_peaks_M,
        n_clusters_sigma,
        n_clusters_M,
        nsigma,
        nM,
        flagged,
        fair_filename,
    )


# ============================================================
# Main execution
# ============================================================

def main():
    """Main entry point for internal MMR analysis with ML-FAIR.

    Loads asteroid orbital elements, runs N-body simulations for each asteroid
    using REBOUND, analyzes resonant angles using KDE and OPTICS clustering,
    and generates summary results comparing predicted and ground-truth resonance
    counts. Outputs include visualization plots and CSV summary files.
    """
    print(f"Working directory: {os.getcwd()}")
    asteroides = ler_elementos_particles()
    print(f"{len(asteroides)} asteroids loaded.")

    real_map = load_real_results("real")
    if real_map:
        print(f"Loaded ground-truth labels for {len(real_map)} asteroids from file: real")
    else:
        print("No ground-truth file found (or it was empty). Metrics will be skipped.")

    args_list = list(enumerate(asteroides, start=0))

    from tqdm import tqdm
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
        results = list(tqdm(executor.map(simular_asteroide, args_list), total=len(args_list)))

    results.sort(key=lambda x: x[0])

    summary_rows = []
    a_planet = float(all_results[PLANET_CFG["all_results_index"]]["a"])
    resonance_column = f"a_res_{SELECTED_PLANET}_AU"

    for (
        _, name, a_input, _log_lines, filename,
        n_peaks_s, n_peaks_M, n_clusters_sigma, n_clusters_M,
        nsigma, nM, flagged, fair_filename,
    ) in results:
        if flagged:
            print("\n" + "=" * 60)
            print(f"FLAGGED asteroid: {name}")
            print(f"Figure: {filename}")
            if fair_filename is not None:
                print(f"FAIR-only figure: {fair_filename}")
                print(
                    f"n_peaks_s={n_peaks_s}, n_peaks_M={n_peaks_M}, "
                    f"n_clusters_sigma={n_clusters_sigma}, n_clusters_M={n_clusters_M}"
                )

                fig_handle = None
                if fair_filename is not None and os.path.exists(fair_filename):
                    fig_handle = show_image_matplotlib(fair_filename)
                elif os.path.exists(filename):
                    fig_handle = show_image_matplotlib(filename)

                while True:
                    try:
                        nsigma = int(input("Enter manual nsigma (integer): ").strip())
                        nM = int(input("Enter manual nM (integer): ").strip())
                        break
                    except ValueError:
                        print("Invalid input. Please enter integers (e.g., 0, 1, 2...).")

                if fig_handle is not None:
                    plt.close(fig_handle)

        # Inner-resonance FAIR mapping:
        #   vertical intersections = p + q  (here: nsigma)
        #   horizontal intersections = q    (here: nM)
        #   => p = (p+q) - q
        p_plus_q = nsigma
        q = nM

        p = np.nan
        pq = np.nan
        Delta_a = np.nan
        a_res_au = np.nan

        if np.isfinite(p_plus_q) and np.isfinite(q):
            p_plus_q = int(p_plus_q)
            q = int(q)
            p = p_plus_q - q
            pq = p_plus_q

        if np.isfinite(p) and np.isfinite(pq) and (p > 0) and (q >= 0):
            a_res_au = resonance_a_inner(a_planet, int(p), int(pq))

            if np.isfinite(a_res_au):
                Delta_a = abs(a_res_au - a_input)
                if Delta_a < 0.4:
                    summary_rows.append({
                        "asteroid": name,
                        "planet": SELECTED_PLANET,
                        "p": int(p),
                        "pq": int(pq),
                        "a_input_AU": f"{a_input:.4g}",
                        resonance_column: f"{a_res_au:.4g}",
                        "Delta_a": f"{Delta_a:.4f}",
                        "flagged": bool(flagged),
                    })

    if summary_rows:
        df_summary = pd.DataFrame(summary_rows)
        summary_file = f"resonant_summary_internal_{SELECTED_PLANET}.csv"
        df_summary.to_csv(summary_file, index=False)
        print(f"\nSaved clustering summary to: {summary_file}")
    else:
        print("\nNo resonant asteroids found (Delta_a < 0.4). No summary file written.")

    print("\n=== Results by asteroid ===")
    for (
        _, name, a_input, log_lines, filename,
        _n_peaks_s, _n_peaks_M, _n_clusters_sigma, _n_clusters_M,
        _nsigma, _nM, flagged, fair_filename,
    ) in results:
        print("\n" + "\n".join(log_lines))
        print(f"Output figure: {filename}")

    if real_map:
        y_true_sigma, y_pred_sigma = [], []
        y_true_M, y_pred_M = [], []
        matched = 0

        for (
            _, name, a_input, _, _,
            _n_peaks_s, _n_peaks_M, n_clusters_sigma, n_clusters_M,
            _nsigma, _nM, flagged, fair_filename,
        ) in results:
            if name in real_map:
                matched += 1
                y_true_sigma.append(real_map[name]["n_clusters_sigma"])
                y_pred_sigma.append(int(n_clusters_sigma))
                y_true_M.append(real_map[name]["n_clusters_M"])
                y_pred_M.append(int(n_clusters_M))

        print("\n=== Evaluation against ground truth ===")
        print(f"Matched {matched} / {len(real_map)} asteroids by name.")

        if matched > 0:
            print(f"σ accuracy: {accuracy_score(y_true_sigma, y_pred_sigma):.3f}")
            print(f"M accuracy: {accuracy_score(y_true_M, y_pred_M):.3f}")

    print("\nSimulations completed.")


if __name__ == "__main__":
    main()
