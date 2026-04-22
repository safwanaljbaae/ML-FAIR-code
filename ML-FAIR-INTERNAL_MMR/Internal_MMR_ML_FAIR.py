"""Internal MMR ML-FAIR - Analysis of internal mean-motion resonances."""
import argparse
import concurrent.futures
import os

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rebound
from angles import normalize
from tqdm import tqdm

from ml_fair import io, analysis, planets


def compute_sigma_elements_internal(planet, asteroid, primary):
    """Compute the internal resonant angle sigma in the invariable plane.

    Sigma = (Omega + omega + M)_planet - (Omega + omega + M)_asteroid
    """
    orb_p = planet.orbit(primary=primary)
    orb_a = asteroid.orbit(primary=primary)

    lam_p = (orb_p.Omega + orb_p.omega + orb_p.M) % (2 * np.pi)
    lam_a = (orb_a.Omega + orb_a.omega + orb_a.M) % (2 * np.pi)

    return (lam_p - lam_a) % (2 * np.pi)


def simular_asteroide(args):
    """Run N-body simulation for a single asteroid and analyze resonant behavior."""
    idx, elementos_asteroide = args
    a_ast, e_ast, inc_ast, Omega_ast, omega_ast, M_ast = elementos_asteroide

    safe_idx = min(idx, len(io.ast_id) - 1)
    name = io.ast_id[safe_idx]
    a_input = io.ast_a_input[safe_idx]

    log_lines = [name]
    flagged = False
    nsigma = np.nan
    nM = np.nan

    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "bs"
    sim.add(m=0.9999999999950272)

    i = 0
    for row in planets.ALL_RESULTS[0:5]:
        sim.add(
            m=planets.MASSAS[i],
            a=row["a"], e=row["e"],
            inc=np.radians(row["incl"]), Omega=np.radians(row["Omega"]),
            omega=np.radians(row["w"]), M=np.radians(row["M"]),
            primary=sim.particles[0],
        )
        i += 1

    sim.add(
        m=0, a=a_ast, e=e_ast, inc=inc_ast, Omega=Omega_ast,
        omega=omega_ast, M=M_ast, primary=sim.particles[0],
    )

    for row in planets.ALL_RESULTS[5:]:
        sim.add(
            m=planets.MASSAS[i],
            a=row["a"], e=row["e"],
            inc=np.radians(row["incl"]), Omega=np.radians(row["Omega"]),
            omega=np.radians(row["w"]), M=np.radians(row["M"]),
            primary=sim.particles[0],
        )
        i += 1

    sim.move_to_com()

    times = np.arange(0.0, 410000 / 365.25 + 0.1095, 0.1095)
    N = len(times)
    sigma_time = np.zeros(N)
    M_Ast = np.zeros(N)

    for k, t in enumerate(times):
        sim.integrate(t)
        orb = sim.particles[6].orbit(primary=sim.particles[0])
        M_Ast[k] = np.degrees(normalize(orb.M))
        sigma_rad = compute_sigma_elements_internal(
            sim.particles[planets.PLANET_CFG["sim_index"]],
            sim.particles[6], sim.particles[0],
        )
        sigma_time[k] = np.degrees(normalize(sigma_rad))

    M = np.mod(M_Ast, 360.0)
    sigma = np.mod(sigma_time, 360.0)

    tol, dens_cov = 5.0, 0.945
    mask_M0 = (M <= tol) | (M >= 360.0 - tol)
    sigma_at_M0_list = sigma[mask_M0]
    mask_sigma0 = (sigma <= tol) | (sigma >= 360.0 - tol)
    M_at_sigma0_list = M[mask_sigma0]

    frac_sigma, _, _ = analysis.circular_union_coverage(sigma_at_M0_list, halfwidth_deg=2)
    frac_M, _, _ = analysis.circular_union_coverage(M_at_sigma0_list, halfwidth_deg=2)

    log_lines.append(f"σ@M≈0 coverage fraction: {frac_sigma}")
    log_lines.append(f"M@σ≈0 coverage fraction: {frac_M}")

    if (frac_sigma >= dens_cov) and (frac_M >= dens_cov):
        n_peaks_s = n_peaks_M = 0
        n_clusters_sigma = n_clusters_M = 0
        peak_pos_s = peak_pos_M = np.array([])
        grid_s = grid_M = np.array([])
        kde_vals_s = kde_vals_M = np.array([])
        log_lines.append(f"High density coverage (>= {dens_cov}); skipping KDE and OPTICS.")
    else:
        n_peaks_s, peak_pos_s, grid_s, kde_vals_s = analysis.count_peaks_circular(
            sigma_at_M0_list, bandwidth_deg=2.5, prom_frac_range=0.40,
        )
        log_lines.append(f"Number of peaks for σ at M≈0: {n_peaks_s}")

        n_peaks_M, peak_pos_M, grid_M, kde_vals_M = analysis.count_peaks_circular(
            M_at_sigma0_list, bandwidth_deg=5, prom_frac_range=0.40,
        )
        log_lines.append(f"Number of peaks for M at σ≈0: {n_peaks_M}")

        _, n_clusters_sigma = analysis.optics_circular_angles(sigma_at_M0_list, max_eps_deg=15, min_samples=2)
        log_lines.append(f"Number of clusters in σ at M≈0: {n_clusters_sigma}")

        _, n_clusters_M = analysis.optics_circular_angles(M_at_sigma0_list, max_eps_deg=15, min_samples=2)
        log_lines.append(f"Number of clusters in M at σ≈0: {n_clusters_M}")

        if n_peaks_s == 1 and n_peaks_M == 1 and n_clusters_sigma == 1 and n_clusters_M == 1:
            n_peaks_s = n_peaks_M = 0
            n_clusters_sigma = n_clusters_M = 0
            nsigma = nM = 0
            log_lines.append("Single peak/cluster in both → resetting to 0.")

        if (n_peaks_s > 0 and n_clusters_sigma > 0 and n_peaks_s != n_clusters_sigma) or \
           (n_peaks_M > 0 and n_clusters_M > 0 and n_peaks_M != n_clusters_M):
            flagged = True
            log_lines.append("FLAGGED: KDE peaks and OPTICS clusters disagree.")

    if (n_peaks_s == n_clusters_sigma) and (n_peaks_M == n_clusters_M):
        nsigma = int(n_peaks_s)
        nM = int(n_peaks_M)
        log_lines.append(f"Agreed counts → nsigma={nsigma}, nM={nM}")
    else:
        nsigma = np.nan
        nM = np.nan

    cvals = np.linspace(0.0, 1.0, len(M_Ast))
    fair_filename = None

    if flagged:
        fair_filename = f"res_forgacs_{planets.SELECTED_PLANET}_{name}_FAIR_FLAGGED.png"
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

    axs[0].scatter(M_Ast, sigma_time, c=cvals, cmap="coolwarm", marker=".", s=1)
    axs[0].set_xlabel("M_ast [deg]", fontsize=16)
    axs[0].set_ylabel(r"$\sigma$ [deg]", fontsize=14)
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
        axs[3].scatter(peak_pos_s, kde_vals_s[np.searchsorted(grid_s, peak_pos_s)], color="r", zorder=3)
    axs[3].set_xlabel(r"$\sigma$ at $M \approx 0^\circ$ [deg]")
    axs[3].set_ylabel("KDE")

    axs[4].plot(grid_M, kde_vals_M, color="k")
    if peak_pos_M.size:
        axs[4].scatter(peak_pos_M, kde_vals_M[np.searchsorted(grid_M, peak_pos_M)], color="r", zorder=3)
    axs[4].set_xlabel(r"$M$ at $\sigma \approx 0^\circ$ [deg]")
    axs[4].set_ylabel("KDE")

    plt.tight_layout()
    filename = f"res_forgacs_{planets.SELECTED_PLANET}_{name}{'_FLAGGED' if flagged else ''}.png"
    plt.savefig(filename)
    plt.close()

    return (idx, name, a_input, log_lines, filename, n_peaks_s, n_peaks_M,
            n_clusters_sigma, n_clusters_M, nsigma, nM, flagged, fair_filename)


def main():
    """Main entry point for internal MMR analysis."""
    print(f"Working directory: {os.getcwd()}")

    planets.ALL_RESULTS = planets.get_planets()

    io.ast_id, io.ast_a_input = io.load_asteroid_list("ast_VRI.txt")
    asteroides = io.ler_elementos_particles()
    print(f"{len(asteroides)} asteroids loaded.")

    real_map = io.load_real_results("real")
    if real_map:
        print(f"Loaded ground-truth labels for {len(real_map)} asteroids.")
    else:
        print("No ground-truth file found. Metrics will be skipped.")

    parser = argparse.ArgumentParser(description="Compute internal MMR locations.")
    parser.add_argument("--planet", type=str, help="Planet name")
    args = parser.parse_args()

    if args.planet:
        planet_name = args.planet.strip().lower()
    else:
        valid = ", ".join(planets.PLANET_INFO.keys())
        planet_name = input(f"Choose a planet ({valid}): ").strip().lower()

    planets.PLANET_CFG = planets.get_simulation_config(planet_name)
    planets.SELECTED_PLANET = planet_name
    print(f"Selected planet: {planets.PLANET_CFG['label']}")

    args_list = list(enumerate(asteroides))

    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
        results = list(tqdm(executor.map(simular_asteroide, args_list), total=len(args_list)))

    results.sort(key=lambda x: x[0])

    summary_rows = []
    a_planet = float(planets.ALL_RESULTS[planets.PLANET_CFG["all_results_index"]]["a"])
    resonance_column = f"a_res_{planets.SELECTED_PLANET}_AU"

    for (idx, name, a_input, _log, _fname, _nps, _npm, ncs, ncm, ns, nM, flagged, _ffn) in results:
        if flagged:
            print(f"\n=== FLAGGED: {name} ===")
            nsigma = int(input("Enter nsigma (integer): ").strip())
            nM = int(input("Enter nM (integer): ").strip())

        p, q = int(ns) if np.isfinite(ns) else 0, int(nM) if np.isfinite(nM) else 0
        if p > 0 and q >= 0:
            pq = p + q
            a_res_au = analysis.resonance_a_inner(a_planet, p, pq)
            if np.isfinite(a_res_au):
                Delta_a = abs(a_res_au - a_input)
                if Delta_a < 0.4:
                    summary_rows.append({
                        "asteroid": name, "planet": planets.SELECTED_PLANET,
                        "p": p, "pq": pq, "a_input_AU": f"{a_input:.4g}",
                        resonance_column: f"{a_res_au:.4g}", "Delta_a": f"{Delta_a:.4f}",
                        "flagged": bool(flagged),
                    })

    if summary_rows:
        df = pd.DataFrame(summary_rows)
        df.to_csv(f"resonant_summary_{planets.SELECTED_PLANET}.csv", index=False)
        print(f"\nSaved summary to resonant_summary_{planets.SELECTED_PLANET}.csv")

    if real_map:
        y_true, y_pred = [], []
        for (_, name, _, _, _, _, _, ncs, ncm, _, _, _, _) in results:
            if name in real_map:
                y_true.append(real_map[name]["n_clusters_sigma"])
                y_pred.append(int(ncs))
        from sklearn.metrics import accuracy_score
        print(f"\nAccuracy: {accuracy_score(y_true, y_pred):.3f}")

    print("\nSimulations completed.")


if __name__ == "__main__":
    main()