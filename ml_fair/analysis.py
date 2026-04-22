"""Analysis utilities for ML-FAIR resonance detection."""
import numpy as np
from scipy.ndimage import maximum_filter1d
from scipy.signal import peak_prominences
from scipy.stats import gaussian_kde
from sklearn.cluster import OPTICS


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


def resonance_a_outer(a_planet_au: float, p: int, pq: int) -> float:
    """Calculate the semi-major axis of an outer p:q mean motion resonance.

    Uses Kepler's third law to compute the resonance location:
    a_res = a_planet * (pq/p)^(2/3).

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
    return a_planet_au * (pq / p) ** (2.0 / 3.0)


def resonance_a_inner(a_planet_au: float, p: int, pq: int) -> float:
    """Calculate the semi-major axis of an inner p:q mean motion resonance.

    Uses Kepler's third law to compute the resonance location:
    a_res = a_planet * (p/pq)^(2/3).

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