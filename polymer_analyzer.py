import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os
import json


def get_param_from_state(wm, category, key):
    """Helper to safely read a parameter from state.json."""
    try:
        with open(wm.state_file, 'r') as f:
            state = json.load(f)
        return state.get(category, {}).get(key)
    except Exception:
        return None


def _safe_time_array(traj):
    """Return a usable time array for plotting and metadata."""
    if getattr(traj, "time", None) is not None and len(traj.time) == traj.n_frames:
        return np.asarray(traj.time, dtype=float)
    return np.arange(traj.n_frames, dtype=float)


def _is_hydrogen(atom):
    """Robust hydrogen detection for topologies where elements may be missing."""
    try:
        if atom.element is not None and atom.element.symbol == "H":
            return True
    except Exception:
        pass
    name = getattr(atom, "name", "") or ""
    return name.strip().upper().startswith("H")


def _select_end_to_end_pair(traj, use_heavy_atoms=True):
    """
    Select terminal atoms for an end-to-end distance estimate.

    The previous implementation used atom 0 and atom n_atoms-1. For generated
    polymer structures that include hydrogens, this may select terminal hydrogen
    atoms rather than the terminal heavy atoms of the chain. As a safer default,
    we use the first and last non-hydrogen atoms in the topology order, falling
    back to the first and last atom only if heavy atoms cannot be identified.
    """
    if use_heavy_atoms:
        heavy_indices = [atom.index for atom in traj.topology.atoms if not _is_hydrogen(atom)]
        if len(heavy_indices) >= 2:
            return int(heavy_indices[0]), int(heavy_indices[-1]), "first_last_heavy_atoms"

    return 0, int(traj.n_atoms - 1), "first_last_atoms_fallback"


def _analysis_start_frame(n_frames, equilibration_fraction):
    """Convert an equilibration fraction into a valid trajectory start frame."""
    if n_frames <= 0:
        raise ValueError("Trajectory contains no frames.")

    equilibration_fraction = float(equilibration_fraction)
    if not 0.0 <= equilibration_fraction < 1.0:
        raise ValueError("equilibration_fraction must be >= 0.0 and < 1.0")

    start_frame = int(np.floor(n_frames * equilibration_fraction))

    # Ensure at least one frame remains for very short trajectories.
    if start_frame >= n_frames:
        start_frame = n_frames - 1

    return start_frame


def analyze_trajectory(wm, equilibration_fraction=0.5, use_heavy_atom_endpoints=True):
    """
    Calculates Rg and End-to-End distance with error analysis.

    By default, statistics are computed only on the second half of the trajectory
    (equilibration_fraction=0.5), while the full trajectory is still plotted.
    This avoids mixing the initial non-equilibrated relaxation frames with the
    demonstrative production-window estimate.

    Results are saved to results.json and a Matplotlib figure is returned.
    """
    # 1. Retrieve paths & metadata
    traj_path = wm.get_path("trajectory")
    top_path = wm.get_path("gro_file")

    n_monomers = get_param_from_state(wm, "parameters", "target_degree_polymerization")

    if not traj_path or not os.path.exists(traj_path):
        print("❌ Error: Trajectory file not found.")
        return None
    if not top_path or not os.path.exists(top_path):
        print("❌ Error: Topology/GRO file not found.")
        return None

    # 2. Load Trajectory
    print(f"   -> Loading trajectory: {os.path.basename(traj_path)}...")
    traj = md.load(traj_path, top=top_path)

    if traj.n_frames < 2:
        raise ValueError(
            "The trajectory contains fewer than 2 saved frames. "
            "Increase the MD length or decrease the reporting interval."
        )

    time_ps = _safe_time_array(traj)
    start_frame = _analysis_start_frame(traj.n_frames, equilibration_fraction)
    analysis_slice = slice(start_frame, None)
    frames_used = traj.n_frames - start_frame

    # ---------------------------------------------------------
    # CALCULATION 1: Radius of Gyration (Rg)
    # ---------------------------------------------------------
    rg_all = md.compute_rg(traj)
    rg = rg_all[analysis_slice]
    avg_rg = np.mean(rg)
    std_rg = np.std(rg)

    # ---------------------------------------------------------
    # CALCULATION 2: End-to-End Distance (Ree)
    # ---------------------------------------------------------
    idx_start, idx_end, endpoint_method = _select_end_to_end_pair(
        traj, use_heavy_atoms=use_heavy_atom_endpoints
    )
    pairs = np.array([[idx_start, idx_end]])
    ree_all = md.compute_distances(traj, pairs)[:, 0]
    ree = ree_all[analysis_slice]
    avg_ree = np.mean(ree)
    std_ree = np.std(ree)

    analysis_note = (
        f"Statistics computed on the last {frames_used}/{traj.n_frames} saved frames "
        f"after excluding the first {equilibration_fraction:.0%} of the trajectory."
    )

    # ---------------------------------------------------------
    # SAVE RESULTS TO JSON
    # ---------------------------------------------------------
    wm.save_result("analysis_window", {
        "description": "second_half_of_trajectory" if equilibration_fraction == 0.5 else "post_equilibration_window",
        "equilibration_fraction_excluded": float(equilibration_fraction),
        "start_frame_index": int(start_frame),
        "frames_used": int(frames_used),
        "total_frames": int(traj.n_frames),
        "start_time_ps": round(float(time_ps[start_frame]), 4),
        "end_time_ps": round(float(time_ps[-1]), 4),
        "note": analysis_note,
    })

    wm.save_result("radius_of_gyration", {
        "mean": round(float(avg_rg), 4),
        "std": round(float(std_rg), 4),
        "units": "nm",
        "analysis_window": "second half of saved trajectory" if equilibration_fraction == 0.5 else "post-equilibration window",
        "note": "Mean and standard deviation exclude the initial non-equilibrium part of the trajectory.",
    })

    wm.save_result("end_to_end_distance", {
        "mean": round(float(avg_ree), 4),
        "std": round(float(std_ree), 4),
        "units": "nm",
        "analysis_window": "second half of saved trajectory" if equilibration_fraction == 0.5 else "post-equilibration window",
        "endpoint_method": endpoint_method,
        "endpoint_atom_indices": [int(idx_start), int(idx_end)],
        "note": "Mean and standard deviation exclude the initial non-equilibrium part of the trajectory.",
    })

    if n_monomers:
        wm.save_result("degree_of_polymerization", int(n_monomers))

    print("✅ Calculations Complete.")
    print(f"   Monomers: {n_monomers}")
    print(f"   Analysis window: frames {start_frame}–{traj.n_frames - 1} "
          f"({frames_used}/{traj.n_frames} frames)")
    print(f"   Time window: {time_ps[start_frame]:.3f}–{time_ps[-1]:.3f} ps")
    print(f"   End-to-end atoms: {idx_start} → {idx_end} ({endpoint_method})")
    print(f"   Rg:  {avg_rg:.3f} ± {std_rg:.3f} nm")
    print(f"   Ree: {avg_ree:.3f} ± {std_ree:.3f} nm")

    # ---------------------------------------------------------
    # PLOTTING
    # ---------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    # Helper for excluded region shading.
    def shade_excluded(axis):
        if start_frame > 0:
            axis.axvspan(time_ps[0], time_ps[start_frame], alpha=0.12,
                         label="Excluded initial window")

    # --- Plot Rg ---
    ax[0].plot(time_ps, rg_all, color='#00796B', linewidth=1.5, label='$R_g$')
    shade_excluded(ax[0])
    ax[0].fill_between(time_ps[start_frame:], avg_rg - std_rg, avg_rg + std_rg,
                       color='#00796B', alpha=0.2,
                       label='Production mean ± std')
    ax[0].hlines(avg_rg, time_ps[start_frame], time_ps[-1],
                 color='black', linestyle='--', alpha=0.8,
                 label=f'Production mean: {avg_rg:.2f}')
    ax[0].set_title("Compactness ($R_g$)", fontsize=14)
    ax[0].set_xlabel("Time (ps)")
    ax[0].set_ylabel("Radius (nm)")
    ax[0].legend(loc='best')
    ax[0].grid(True, alpha=0.3)

    # --- Plot End-to-End ---
    ax[1].plot(time_ps, ree_all, color='#512DA8', linewidth=1.5, label='$R_{ee}$')
    shade_excluded(ax[1])
    ax[1].fill_between(time_ps[start_frame:], avg_ree - std_ree, avg_ree + std_ree,
                       color='#512DA8', alpha=0.2,
                       label=r'production $\bar{x}\pm\sigma$')
    ax[1].hlines(avg_ree, time_ps[start_frame], time_ps[-1],
                 color='black', linestyle='--', alpha=0.8,
                 label=f'Production mean: {avg_ree:.2f}')
    ax[1].set_title("Extension ($R_{ee}$)", fontsize=14)
    ax[1].set_xlabel("Time (ps)")
    ax[1].set_ylabel("Distance (nm)")
    ax[1].legend(loc='best')
    ax[1].grid(True, alpha=0.3)

    fig.suptitle("Structural descriptors: statistics computed after initial equilibration window", y=1.03)
    plt.tight_layout()
    return fig
