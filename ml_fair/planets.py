"""Planet data and configuration for ML-FAIR simulations."""
import numpy as np
import astropy.table
import astropy.units as u
from astroquery.jplhorizons import Horizons


def _to_invariant_plane(a, e, inc, Omega, omega):
    """Transform orbital elements from ecliptic of J2000 to invariable plane."""
    try:
        from astropy.coordinates.matrix_utilities import rotation_matrix
    except ImportError:
        from astropy.coordinates.angles import rotation_matrix

    I_INV = 1.5833333 * u.deg
    OMEGA_INV = 107.0 * u.deg

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

MASSAS = [
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

ALL_RESULTS = None


def get_planets(force_refresh=False):
    """Get planetary orbital elements from JPL Horizons.

    Parameters
    ----------
    force_refresh : bool, optional
        Force re-fetching from Horizons even if cached. Default is False.

    Returns
    -------
    astropy.table.Table
        Table of planetary orbital elements transformed to invariable plane.
    """
    global ALL_RESULTS
    if ALL_RESULTS is not None and not force_refresh:
        return ALL_RESULTS

    bodies = ["199", "299", "399", "301", "499", "599", "699", "799", "899"]
    all_results = astropy.table.vstack([
        Horizons(id=i, location="500@10", epochs=2460800.5).elements()
        for i in bodies
    ])

    for row in all_results:
        a, e = row["a"], row["e"]
        inc, Omega, omega = row["incl"], row["Omega"], row["w"]
        a, e, inc_new, Omega_new, omega_new = _to_invariant_plane(a, e, inc, Omega, omega)
        row["incl"] = inc_new
        row["Omega"] = Omega_new
        row["w"] = omega_new

    ALL_RESULTS = all_results
    return all_results


def get_simulation_config(planet_name):
    """Get simulation configuration for a specific planet.

    Parameters
    ----------
    planet_name : str
        Name of the planet (e.g., 'mercury', 'earth').

    Returns
    -------
    dict
        Configuration dictionary with sim_index, all_results_index, and label.

    Raises
    ------
    ValueError
        If planet_name is not valid.
    """
    planet_name = planet_name.lower()
    if planet_name not in PLANET_INFO:
        raise ValueError(
            f"Invalid planet '{planet_name}'. Choose one of: {', '.join(PLANET_INFO.keys())}"
        )
    return PLANET_INFO[planet_name]