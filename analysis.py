import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os

def analyze_trajectory(wm):
    """
    Calculates Rg and End-to-End distance from the trajectory stored in the WorkflowManager.
    Saves results to results.json and returns a Matplotlib figure.
    
    Args:
        wm (WorkflowManager): The active workflow manager instance containing paths.
    
    Returns:
        fig (matplotlib.figure.Figure): The plot object to display in the notebook.
    """
    # 1. Retrieve paths from State
    traj_path = wm.get_path("trajectory")
    top_path = wm.get_path("gro_file")

    if not traj_path or not os.path.exists(traj_path):
        print("❌ Error: Trajectory file not found. Did the simulation run?")
        return None

    # 2. Load Trajectory
    # We use the GRO file as topology because it contains atom names
    print(f"   -> Loading trajectory: {os.path.basename(traj_path)}...")
    traj = md.load(traj_path, top=top_path)
    
    # ---------------------------------------------------------
    # CALCULATION 1: Radius of Gyration (Rg)
    # Measures how compact the polymer is (low = globular, high = extended)
    # ---------------------------------------------------------
    rg = md.compute_rg(traj)
    avg_rg = np.mean(rg)
    
    # ---------------------------------------------------------
    # CALCULATION 2: End-to-End Distance (Ree)
    # Measures the distance between the first and last atom
    # ---------------------------------------------------------
    first_atom = 0
    last_atom = traj.n_atoms - 1
    pairs = np.array([[first_atom, last_atom]])
    
    # compute_distances returns shape (n_frames, n_pairs), we want 1D array
    ree = md.compute_distances(traj, pairs)[:, 0] 
    avg_ree = np.mean(ree)

    # ---------------------------------------------------------
    # SAVE RESULTS TO JSON
    # ---------------------------------------------------------
    wm.save_result("radius_of_gyration_nm", round(float(avg_rg), 4), units="nm")
    wm.save_result("end_to_end_distance_nm", round(float(avg_ree), 4), units="nm")

    print(f"✅ Calculations Complete.")
    print(f"   Avg Radius of Gyration:  {avg_rg:.3f} nm")
    print(f"   Avg End-to-End Distance: {avg_ree:.3f} nm")

    # ---------------------------------------------------------
    # PLOTTING
    # ---------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    # Plot Rg
    ax[0].plot(traj.time, rg, color='teal', linewidth=1.5, label='$R_g$')
    ax[0].axhline(avg_rg, color='black', linestyle='--', alpha=0.7, label=f'Mean: {avg_rg:.2f} nm')
    ax[0].set_title("Compactness ($R_g$)", fontsize=14)
    ax[0].set_xlabel("Time (ps)")
    ax[0].set_ylabel("Radius (nm)")
    ax[0].legend()
    ax[0].grid(True, alpha=0.3)

    # Plot End-to-End
    ax[1].plot(traj.time, ree, color='purple', linewidth=1.5, label='$R_{ee}$')
    ax[1].axhline(avg_ree, color='black', linestyle='--', alpha=0.7, label=f'Mean: {avg_ree:.2f} nm')
    ax[1].set_title("Extension ($R_{ee}$)", fontsize=14)
    ax[1].set_xlabel("Time (ps)")
    ax[1].set_ylabel("Distance (nm)")
    ax[1].legend()
    ax[1].grid(True, alpha=0.3)

    plt.tight_layout()
    
    return fig
