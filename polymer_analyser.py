import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import io
import contextlib


# =============================================================================
# Helper functions
# =============================================================================


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


def _safe_correlation(x, y):
    """
    Return the Pearson correlation between two arrays.

    If one of the arrays is constant, the correlation is not defined and NaN is
    returned instead of raising a warning/error.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) < 2 or len(y) < 2:
        return float("nan")

    if np.std(x) == 0.0 or np.std(y) == 0.0:
        return float("nan")

    return float(np.corrcoef(x, y)[0, 1])


def _zscore(values):
    """
    Return z-scored values for visual comparison.

    If the signal is constant, return zeros.
    """

    values = np.asarray(values, dtype=float)
    std = np.std(values)

    if std == 0.0:
        return np.zeros_like(values)

    return (values - np.mean(values)) / std


# =============================================================================
# Atom selection and endpoint selection
# =============================================================================


def _select_analysis_atoms(traj, atom_selection=None):
    """
    Select atoms used for polymer analysis.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Loaded trajectory.
    atom_selection : str, array-like or None
        - None: use all atoms, preserving the previous behaviour.
        - str: MDTraj atom selection string, e.g. "resname POL".
        - array-like: explicit atom indices.

    Returns
    -------
    atom_indices : np.ndarray
        Selected atom indices in the original topology.
    selection_description : str
        Human-readable selection description.
    """

    if atom_selection is None:
        atom_indices = np.arange(traj.n_atoms, dtype=int)
        selection_description = "all atoms"

    elif isinstance(atom_selection, str):
        try:
            atom_indices = traj.topology.select(atom_selection)
        except Exception as exc:
            raise ValueError(
                f"Invalid MDTraj atom selection: {atom_selection!r}"
            ) from exc

        selection_description = atom_selection

    else:
        atom_indices = np.asarray(atom_selection, dtype=int)
        selection_description = "explicit atom indices"

    atom_indices = np.asarray(atom_indices, dtype=int)

    if atom_indices.ndim != 1:
        raise ValueError("atom_selection must define a one-dimensional atom index list.")

    if len(atom_indices) < 2:
        raise ValueError(
            f"Atom selection {selection_description!r} contains fewer than 2 atoms."
        )

    if np.any(atom_indices < 0) or np.any(atom_indices >= traj.n_atoms):
        raise ValueError("atom_selection contains atom indices outside the trajectory.")

    # Remove duplicates while preserving order.
    _, unique_positions = np.unique(atom_indices, return_index=True)
    atom_indices = atom_indices[np.sort(unique_positions)]

    return atom_indices, selection_description


def _select_end_to_end_pair(traj, atom_indices=None, use_heavy_atoms=True):
    """
    Select terminal atoms for the end-to-end distance within a selected subset.

    This is the key change relative to the previous version: Ree is no longer
    computed between the first and last atom of the whole system unless the
    selected subset is explicitly the whole system.

    By default, the first and last heavy atoms inside the selected subset are
    used. This is usually more meaningful than using hydrogens as endpoints.
    """

    if atom_indices is None:
        atom_indices = np.arange(traj.n_atoms, dtype=int)
    else:
        atom_indices = np.asarray(atom_indices, dtype=int)

    if len(atom_indices) < 2:
        raise ValueError("Atom selection must contain at least two atoms.")

    selected_atoms = [traj.topology.atom(int(i)) for i in atom_indices]

    if use_heavy_atoms:
        heavy_indices = [
            atom.index
            for atom in selected_atoms
            if not _is_hydrogen(atom)
        ]

        if len(heavy_indices) >= 2:
            return (
                int(heavy_indices[0]),
                int(heavy_indices[-1]),
                "first_last_heavy_atoms_in_selection",
            )

    return (
        int(atom_indices[0]),
        int(atom_indices[-1]),
        "first_last_atoms_in_selection_fallback",
    )


def _endpoint_description(endpoint_method, idx_start, idx_end):
    """
    Human-readable description of the atoms used for Ree.
    """

    if endpoint_method == "first_last_heavy_atoms_in_selection":
        return (
            "first-to-last heavy atom distance within the selected atoms "
            f"(atoms {idx_start} → {idx_end})"
        )

    if endpoint_method == "first_last_atoms_in_selection_fallback":
        return (
            "first-to-last atom distance within the selected atoms "
            f"(fallback; atoms {idx_start} → {idx_end})"
        )

    # Kept for backwards-compatible descriptions if older result files are read.
    if endpoint_method == "first_last_heavy_atoms":
        return (
            "first-to-last heavy atom distance "
            f"(atoms {idx_start} → {idx_end})"
        )

    return (
        "first-to-last atom distance "
        f"(fallback; atoms {idx_start} → {idx_end})"
    )


# =============================================================================
# Main analysis function
# =============================================================================


def analyze_trajectory(
    wm,
    equilibration_fraction=0.5,
    atom_selection=None,
    use_heavy_atom_endpoints=True,
    periodic=False,
    trajectory_key="trajectory",
    topology_key="gro_file",
    show_normalized_comparison=True,
    verbose=False,
):
    """
    Analyze polymer trajectory properties.

    Calculates:
    - Radius of gyration, Rg
    - End-to-end distance, Ree
    - Ree/Rg ratio
    - Pearson correlation between Rg and Ree in the analysis window

    Important
    ---------
    Rg and Ree are meaningful polymer conformational descriptors only if the
    selected atoms correspond to one polymer chain. For a bulk system containing
    many chains, solvent, surfaces, nanoparticles, or other components, pass an
    explicit atom_selection for the chain of interest, or compute these
    descriptors chain-by-chain.

    Parameters
    ----------
    wm : workflow manager
        Object exposing get_path() and save_result().
    equilibration_fraction : float, default=0.5
        Fraction of initial frames excluded from summary statistics.
    atom_selection : str, array-like or None, default=None
        Atoms used for Rg and Ree.
        Examples:
        - None: all atoms, backward-compatible but not recommended for bulk.
        - "resname POL": all atoms with residue name POL.
        - "chainid 0 and not element H": chain 0 heavy atoms, if chain IDs exist.
        - [0, 1, 2, ...]: explicit atom indices.
    use_heavy_atom_endpoints : bool, default=True
        Use first and last heavy atoms in the selected atom subset as endpoints.
    periodic : bool, default=False
        Passed to md.compute_distances for Ree.
        Use False if the molecule is already whole/unwrapped.
        Use True only if periodic minimum-image distances are desired.
    trajectory_key : str, default="trajectory"
        Workflow path key for the trajectory.
    topology_key : str, default="gro_file"
        Workflow path key for the topology/GRO file.
    show_normalized_comparison : bool, default=True
        Add a third panel comparing z-scored Rg and Ree trends.
    verbose : bool, default=False
        Print additional details.

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        Matplotlib figure with the analysis plots.
    """

    # ---------------------------------------------------------
    # Retrieve paths and metadata
    # ---------------------------------------------------------

    traj_path = wm.get_path(trajectory_key)
    top_path = wm.get_path(topology_key)

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
        print(f"   -> Loading topology:   {os.path.basename(top_path)}...")

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
    # Select atoms used for polymer descriptors
    # ---------------------------------------------------------

    atom_indices, atom_selection_text = _select_analysis_atoms(
        traj,
        atom_selection=atom_selection,
    )

    if verbose:
        print(f"   -> Atom selection: {atom_selection_text}")
        print(f"   -> Selected atoms: {len(atom_indices)}/{traj.n_atoms}")

    if atom_selection is None and traj.n_atoms > len(atom_indices):
        # This branch is kept for completeness; with atom_selection=None the two
        # values are equal. The explicit warning below is more useful.
        pass

    if atom_selection is None:
        print(
            "⚠️  Atom selection: all atoms. "
            "For bulk/surface/solvated systems, pass atom_selection for one polymer chain."
        )

    # Create a trajectory containing only the selected atoms for Rg.
    # This prevents Rg from being computed over solvent, surfaces, nanoparticles,
    # or multiple chains unless they were explicitly selected.
    traj_selected = traj.atom_slice(atom_indices)

    # ---------------------------------------------------------
    # Radius of gyration
    # ---------------------------------------------------------

    rg_all = md.compute_rg(traj_selected)
    rg = rg_all[analysis_slice]

    avg_rg, std_rg = _mean_and_std(rg)

    # ---------------------------------------------------------
    # End-to-end distance
    # ---------------------------------------------------------

    idx_start, idx_end, endpoint_method = _select_end_to_end_pair(
        traj,
        atom_indices=atom_indices,
        use_heavy_atoms=use_heavy_atom_endpoints,
    )

    endpoint_text = _endpoint_description(
        endpoint_method,
        idx_start,
        idx_end,
    )

    pairs = np.array([[idx_start, idx_end]], dtype=int)
    ree_all = md.compute_distances(traj, pairs, periodic=periodic)[:, 0]
    ree = ree_all[analysis_slice]

    avg_ree, std_ree = _mean_and_std(ree)

    # ---------------------------------------------------------
    # Diagnostics
    # ---------------------------------------------------------

    ree_rg_ratio = float(avg_ree / avg_rg) if avg_rg > 0 else float("nan")
    rg_ree_correlation = _safe_correlation(rg, ree)

    if np.isfinite(rg_ree_correlation) and abs(rg_ree_correlation) > 0.95:
        correlation_note = (
            "Rg and Ree are very strongly correlated in the analysis window. "
            "Check that the atom selection corresponds to one polymer chain and "
            "that the trajectory is not dominated by global/PBC motion."
        )
    else:
        correlation_note = (
            "Correlation computed over the analysis window. High correlation can "
            "be physical for short trajectories, but can also indicate an atom "
            "selection or PBC issue."
        )

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
        "analysis_atom_selection",
        {
            "selection": atom_selection_text,
            "selected_atom_count": int(len(atom_indices)),
            "total_atom_count": int(traj.n_atoms),
            "first_selected_atom_index": int(atom_indices[0]),
            "last_selected_atom_index": int(atom_indices[-1]),
            "note": (
                "Rg and Ree were computed only on this atom selection. "
                "For conformational descriptors, this should usually be one polymer chain."
            ),
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
            "atom_selection": atom_selection_text,
            "note": (
                "Mean and standard deviation exclude the initial "
                "non-equilibrium part of the trajectory. Rg was computed "
                "on the selected atoms only."
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
            "atom_selection": atom_selection_text,
            "endpoint_description": endpoint_text,
            "endpoint_method": endpoint_method,
            "endpoint_atom_indices": [int(idx_start), int(idx_end)],
            "periodic": bool(periodic),
            "note": (
                "Mean and standard deviation exclude the initial "
                "non-equilibrium part of the trajectory. Ree was computed "
                "between terminal atoms within the selected atoms only."
            ),
        },
    )

    _save_result_silent(
        wm,
        "chain_shape_diagnostics",
        {
            "ree_over_rg": round(ree_rg_ratio, 4),
            "rg_ree_correlation": (
                round(rg_ree_correlation, 4)
                if np.isfinite(rg_ree_correlation)
                else None
            ),
            "correlation_note": correlation_note,
            "note": (
                "For an ideal freely-jointed-like chain, Ree/Rg is often around "
                "sqrt(6) ≈ 2.45. Values depend on chain stiffness, length, chemistry, "
                "environment, and sampling."
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

    print("--> Calculations Complete.")
    print(f"   Monomers: {monomer_text}")
    print(f"   Atom selection: {atom_selection_text}")
    print(f"   Selected atoms: {len(atom_indices)}/{traj.n_atoms}")
    print(
        f"   Analysis window: frames {start_frame}–{traj.n_frames - 1} "
        f"({frames_used}/{traj.n_frames} frames)"
    )
    print(
        f"   Time window: {time_ps[start_frame]:.3f}–{time_ps[-1]:.3f} ps"
    )
    print(f"   End-to-end distance: {endpoint_text}")
    print(f"   Rg:        {avg_rg:.3f} ± {std_rg:.3f} nm")
    print(f"   Ree:       {avg_ree:.3f} ± {std_ree:.3f} nm")
    print(f"   Ree/Rg:    {ree_rg_ratio:.3f}")

    if np.isfinite(rg_ree_correlation):
        print(f"   corr(Rg, Ree): {rg_ree_correlation:.3f}")
    else:
        print("   corr(Rg, Ree): not defined")

    if np.isfinite(rg_ree_correlation) and abs(rg_ree_correlation) > 0.95:
        print("   ⚠️  Rg and Ree are very strongly correlated. Check atom selection/PBC.")

    # ---------------------------------------------------------
    # Plotting
    # ---------------------------------------------------------

    if show_normalized_comparison:
        fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    else:
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

    # Normalized comparison plot
    if show_normalized_comparison:
        rg_z = _zscore(rg_all)
        ree_z = _zscore(ree_all)

        ax[2].plot(
            time_ps,
            rg_z,
            color="#00796B",
            linewidth=1.5,
            label="z-score $R_g$",
        )
        ax[2].plot(
            time_ps,
            ree_z,
            color="#512DA8",
            linewidth=1.5,
            label="z-score $R_{ee}$",
        )

        shade_excluded(ax[2])

        ax[2].set_title("Normalized trend comparison", fontsize=14)
        ax[2].set_xlabel("Time (ps)")
        ax[2].set_ylabel("z-score")
        ax[2].legend(loc="best")
        ax[2].grid(True, alpha=0.3)

    fig.suptitle(
        "Polymer structural descriptors",
        y=1.03,
    )

    plt.tight_layout()

    return fig
