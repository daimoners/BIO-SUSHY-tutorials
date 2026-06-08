import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import io
import contextlib


def get_param_from_state(wm, category, key):
    """
    Safely read a parameter from the workflow state.

    Priority:
    1. in-memory wm.state
    2. state.json on disk
    """

    try:
        state = getattr(wm, "state", None)
        if isinstance(state, dict):
            value = state.get(category, {}).get(key)
            if value is not None:
                return value
    except Exception:
        pass

    try:
        with open(wm.state_file, "r") as f:
            state = json.load(f)
        return state.get(category, {}).get(key)
    except Exception:
        return None


def _save_result_silent(wm, key, value, units=None):
    """
    Save a result without printing workflow-level messages.

    This keeps the notebook output focused on the final analysis summary.
    """

    with contextlib.redirect_stdout(io.StringIO()):
        if units is None:
            wm.save_result(key, value)
        else:
            wm.save_result(key, value, units=units)


def _safe_time_array(traj):
    """
    Return a usable time array for plotting and metadata.

    MDTraj stores trajectory time in ps. If missing, frame indices are used.
    """

    if getattr(traj, "time", None) is not None and len(traj.time) == traj.n_frames:
        return np.asarray(traj.time, dtype=float)

    return np.arange(traj.n_frames, dtype=float)


def _is_hydrogen(atom):
    """
    Robust hydrogen detection for topologies where element information
    may be missing.
    """

    try:
        if atom.element is not None and atom.element.symbol == "H":
            return True
    except Exception:
        pass

    name = getattr(atom, "name", "") or ""
    return name.strip().upper().startswith("H")


def _select_end_to_end_pair(traj, use_heavy_atoms=True):
    """
    Select terminal atoms for the end-to-end distance.

    By default, the first and last heavy atoms are used. This is usually more
    meaningful than using the first and last atom in the topology, because those
    may be hydrogens.
    """

    if traj.n_atoms < 2:
        raise ValueError("Trajectory must contain at least two atoms.")

    if use_heavy_atoms:
        heavy_indices = [
            atom.index
            for atom in traj.topology.atoms
            if not _is_hydrogen(atom)
        ]

        if len(heavy_indices) >= 2:
            return (
                int(heavy_indices[0]),
                int(heavy_indices[-1]),
                "first_last_heavy_atoms",
            )

    return 0, int(traj.n_atoms - 1), "first_last_atoms_fallback"


def _endpoint_description(endpoint_method, idx_start, idx_end):
    """
    Human-readable description of the atoms used for Ree.
    """

    if endpoint_method == "first_last_heavy_atoms":
        return (
            "first-to-last heavy atom distance "
            f"(atoms {idx_start} → {idx_end})"
        )

    return (
        "first-to-last atom distance "
        f"(fallback; atoms {idx_start} → {idx_end})"
    )


def _analysis_start_frame(n_frames, equilibration_fraction):
    """
    Convert an equilibration fraction into a valid trajectory start frame.
    """

    if n_frames <= 0:
        raise ValueError("Trajectory contains no frames.")

    equilibration_fraction = float(equilibration_fraction)

    if not 0.0 <= equilibration_fraction < 1.0:
        raise ValueError("equilibration_fraction must be >= 0.0 and < 1.0")

    start_frame = int(np.floor(n_frames * equilibration_fraction))

    if start_frame >= n_frames:
        start_frame = n_frames - 1

    return start_frame


def _mean_and_std(values):
    """
    Return mean and standard deviation as Python floats.
    """

    values = np.asarray(values, dtype=float)
    return float(np.mean(values)), float(np.std(values))


def analyze_trajectory(
    wm,
    equilibration_fraction=0.5,
    use_heavy_atom_endpoints=True,
    verbose=False,
):
    """
    Analyze polymer trajectory properties.

    Calculates:
    - Radius of gyration, Rg
    - End-to-end distance, Ree

    The full trajectory is plotted, but summary statistics are computed only
    after the initial equilibration window. By default, this means the second
    half of the saved trajectory is used.

    Results are saved to results.json and a Matplotlib figure is returned.
    """

    # ---------------------------------------------------------
    # Retrieve paths and metadata
    # ---------------------------------------------------------

    traj_path = wm.get_path("trajectory")
    top_path = wm.get_path("gro_file")

    n_monomers = get_param_from_state(
        wm,
        "parameters",
        "target_degree_polymerization",
    )

    if not traj_path or not os.path.exists(traj_path):
        print("❌ Error: Trajectory file not found.")
        return None

    if not top_path or not os.path.exists(top_path):
        print("❌ Error: Topology/GRO file not found.")
        return None

    # ---------------------------------------------------------
    # Load trajectory
    # ---------------------------------------------------------

    if verbose:
        print(f"   -> Loading trajectory: {os.path.basename(traj_path)}...")

    traj = md.load(traj_path, top=top_path)

    if traj.n_frames < 2:
        raise ValueError(
            "The trajectory contains fewer than 2 saved frames. "
            "Increase the MD length or decrease the reporting interval."
        )

    time_ps = _safe_time_array(traj)

    start_frame = _analysis_start_frame(
        traj.n_frames,
        equilibration_fraction,
    )

    analysis_slice = slice(start_frame, None)
    frames_used = traj.n_frames - start_frame

    window_label = (
        "second half of saved trajectory"
        if equilibration_fraction == 0.5
        else "post-equilibration window"
    )

    analysis_note = (
        f"Statistics computed on the last {frames_used}/{traj.n_frames} "
        f"saved frames after excluding the first "
        f"{equilibration_fraction:.0%} of the trajectory."
    )

    # ---------------------------------------------------------
    # Radius of gyration
    # ---------------------------------------------------------

    rg_all = md.compute_rg(traj)
    rg = rg_all[analysis_slice]

    avg_rg, std_rg = _mean_and_std(rg)

    # ---------------------------------------------------------
    # End-to-end distance
    # ---------------------------------------------------------

    idx_start, idx_end, endpoint_method = _select_end_to_end_pair(
        traj,
        use_heavy_atoms=use_heavy_atom_endpoints,
    )

    endpoint_text = _endpoint_description(
        endpoint_method,
        idx_start,
        idx_end,
    )

    pairs = np.array([[idx_start, idx_end]])
    ree_all = md.compute_distances(traj, pairs)[:, 0]
    ree = ree_all[analysis_slice]

    avg_ree, std_ree = _mean_and_std(ree)

    # ---------------------------------------------------------
    # Save results
    # ---------------------------------------------------------

    _save_result_silent(
        wm,
        "analysis_window",
        {
            "description": (
                "second_half_of_trajectory"
                if equilibration_fraction == 0.5
                else "post_equilibration_window"
            ),
            "equilibration_fraction_excluded": float(equilibration_fraction),
            "start_frame_index": int(start_frame),
            "end_frame_index": int(traj.n_frames - 1),
            "frames_used": int(frames_used),
            "total_frames": int(traj.n_frames),
            "start_time_ps": round(float(time_ps[start_frame]), 4),
            "end_time_ps": round(float(time_ps[-1]), 4),
            "note": analysis_note,
        },
    )

    _save_result_silent(
        wm,
        "radius_of_gyration",
        {
            "mean": round(avg_rg, 4),
            "std": round(std_rg, 4),
            "units": "nm",
            "analysis_window": window_label,
            "note": (
                "Mean and standard deviation exclude the initial "
                "non-equilibrium part of the trajectory."
            ),
        },
    )

    _save_result_silent(
        wm,
        "end_to_end_distance",
        {
            "mean": round(avg_ree, 4),
            "std": round(std_ree, 4),
            "units": "nm",
            "analysis_window": window_label,
            "endpoint_description": endpoint_text,
            "endpoint_method": endpoint_method,
            "endpoint_atom_indices": [int(idx_start), int(idx_end)],
            "note": (
                "Mean and standard deviation exclude the initial "
                "non-equilibrium part of the trajectory."
            ),
        },
    )

    if n_monomers is not None:
        try:
            n_monomers_saved = int(n_monomers)
        except Exception:
            n_monomers_saved = n_monomers

        _save_result_silent(
            wm,
            "degree_of_polymerization",
            n_monomers_saved,
        )

    # ---------------------------------------------------------
    # Compact notebook output
    # ---------------------------------------------------------

    monomer_text = n_monomers if n_monomers is not None else "not available"

    print("✅ Calculations Complete.")
    print(f"   Monomers: {monomer_text}")
    print(
        f"   Analysis window: frames {start_frame}–{traj.n_frames - 1} "
        f"({frames_used}/{traj.n_frames} frames)"
    )
    print(
        f"   Time window: {time_ps[start_frame]:.3f}–{time_ps[-1]:.3f} ps"
    )
    print(f"   End-to-end distance: {endpoint_text}")
    print(f"   Rg:  {avg_rg:.3f} ± {std_rg:.3f} nm")
    print(f"   Ree: {avg_ree:.3f} ± {std_ree:.3f} nm")

    # ---------------------------------------------------------
    # Plotting
    # ---------------------------------------------------------

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    def shade_excluded(axis):
        if start_frame > 0:
            axis.axvspan(
                time_ps[0],
                time_ps[start_frame],
                alpha=0.12,
                label="Excluded initial window",
            )

    # Radius of gyration plot
    ax[0].plot(
        time_ps,
        rg_all,
        color="#00796B",
        linewidth=1.5,
        label="$R_g$",
    )

    shade_excluded(ax[0])

    ax[0].fill_between(
        time_ps[start_frame:],
        avg_rg - std_rg,
        avg_rg + std_rg,
        color="#00796B",
        alpha=0.2,
        label="Analysis mean ± std",
    )

    ax[0].hlines(
        avg_rg,
        time_ps[start_frame],
        time_ps[-1],
        color="black",
        linestyle="--",
        alpha=0.8,
        label=f"Analysis mean: {avg_rg:.2f} nm",
    )

    ax[0].set_title("Compactness ($R_g$)", fontsize=14)
    ax[0].set_xlabel("Time (ps)")
    ax[0].set_ylabel("Radius of gyration (nm)")
    ax[0].legend(loc="best")
    ax[0].grid(True, alpha=0.3)

    # End-to-end distance plot
    ax[1].plot(
        time_ps,
        ree_all,
        color="#512DA8",
        linewidth=1.5,
        label="$R_{ee}$",
    )

    shade_excluded(ax[1])

    ax[1].fill_between(
        time_ps[start_frame:],
        avg_ree - std_ree,
        avg_ree + std_ree,
        color="#512DA8",
        alpha=0.2,
        label="Analysis mean ± std",
    )

    ax[1].hlines(
        avg_ree,
        time_ps[start_frame],
        time_ps[-1],
        color="black",
        linestyle="--",
        alpha=0.8,
        label=f"Analysis mean: {avg_ree:.2f} nm",
    )

    ax[1].set_title("Extension ($R_{ee}$)", fontsize=14)
    ax[1].set_xlabel("Time (ps)")
    ax[1].set_ylabel("End-to-end distance (nm)")
    ax[1].legend(loc="best")
    ax[1].grid(True, alpha=0.3)

    fig.suptitle(
        "Polymer structural descriptors",
        y=1.03,
    )

    plt.tight_layout()

    return fig
