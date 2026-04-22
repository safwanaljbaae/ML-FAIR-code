"""ML-FAIR: Machine-Learning-Assisted Fair Analysis of Mean Motion Resonances.

A Python framework for identifying and analyzing mean-motion resonances (MMRs)
in the Solar System using N-body simulations and clustering techniques.

Submodules
----------
- io: Input/output utilities for asteroid data and orbital elements
- analysis: Core analysis functions for resonance detection
- planets: Planet data and configuration

Usage
-----
from ml_fair import io, analysis, planets

# Load asteroid data
ast_id, ast_a_input = io.load_asteroid_list("ast_VRI.txt")

# Get planet data
planets.get_planets()

# Analyze resonances
n_peaks, peak_pos, grid, kde = analysis.count_peaks_circular(data)
"""

from ml_fair.io import (
    load_asteroid_list,
    load_real_results,
    ler_elementos_particles,
    to_invariant_plane,
)
from ml_fair.analysis import (
    circular_union_coverage,
    count_peaks_circular,
    optics_circular_angles,
    resonance_a_inner,
    resonance_a_outer,
)
from ml_fair.planets import (
    MASSAS,
    PLANET_INFO,
    get_planets,
    get_simulation_config,
)

__version__ = "1.0.0"
__all__ = [
    "io",
    "analysis",
    "planets",
    "load_asteroid_list",
    "load_real_results",
    "ler_elementos_particles",
    "to_invariant_plane",
    "circular_union_coverage",
    "count_peaks_circular",
    "optics_circular_angles",
    "resonance_a_inner",
    "resonance_a_outer",
    "MASSAS",
    "PLANET_INFO",
    "get_planets",
    "get_simulation_config",
]