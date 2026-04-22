"""Common utilities for ML-FAIR resonance analysis."""
import os
import numpy as np
import astropy.table
import astropy.units as u

I_INV = 1.5833333 * u.deg
OMEGA_INV = 107.0 * u.deg


def to_invariant_plane(a, e, inc, Omega, omega):
    """Transform orbital elements from ecliptic of J2000 to invariable plane.

    Parameters
    ----------
    a : float
        Semi-major axis in AU.
    e : float
        Eccentricity.
    inc : float
        Inclination in degrees.
    Omega : float
        Longitude of ascending node in degrees.
    omega : float
        Argument of perihelion in degrees.

    Returns
    -------
    tuple
        (a, e, inc_new, Omega_new, omega_new) - Transformed elements.
    """
    try:
        from astropy.coordinates.matrix_utilities import rotation_matrix
    except ImportError:
        from astropy.coordinates.angles import rotation_matrix

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


def ler_elementos_particles(arquivo="Particles.el", use_invariant=True):
    """Read asteroid orbital elements from a file and transform to invariable plane.

    Parameters
    ----------
    arquivo : str, optional
        Path to the file containing orbital elements. Default is "Particles.el".
    use_invariant : bool, optional
        Whether to transform to invariable plane. Default is True.

    Returns
    -------
    list of tuple
        List of tuples containing (a, e, inc_rad, Omega_rad, omega_rad, M_rad)
        for each asteroid, with angles in radians.
    """
    with open(arquivo, "r") as f:
        linhas = f.readlines()
    elementos = []
    for linha in linhas:
        partes = linha.strip().split()
        a, e, inc, Omega, omega, M = map(float, partes[:6])
        if use_invariant:
            a, e, inc_new, Omega_new, omega_new = to_invariant_plane(a, e, inc, Omega, omega)
            elementos.append((
                a, e, np.radians(inc_new), np.radians(Omega_new),
                np.radians(omega_new), np.radians(M),
            ))
        else:
            elementos.append((
                a, e, np.radians(inc), np.radians(Omega),
                np.radians(omega), np.radians(M),
            ))
    return elementos


def load_asteroid_list(input_file="ast_VRI.txt"):
    """Load asteroid names and input semi-major axes from a file.

    Parameters
    ----------
    input_file : str, optional
        Path to the asteroid list file.

    Returns
    -------
    tuple
        (ast_id, ast_a_input) - List of asteroid names and corresponding semi-major axes.
    """
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
    return ast_id, ast_a_input