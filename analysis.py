import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import os
import json

def get_param_from_state(wm, category, key):
    """Helper to safely read a parameter from state.json"""
    try:
        with open(wm.state_file, 'r') as f:
            state = json.load(f)
        return state.get(category, {}).get(key)
    except:
        return None

def analyze_trajectory(wm):
    """
    Calculates Rg and End-to-End distance with Error analysis.
    Saves results to results.json and returns a Matplotlib figure.
    """
    # 1. Retrieve paths & metadata
    traj_path = wm.get_path("trajectory")
    top_path = wm.get_path("gro_file")
    
    # Retrieve monomer count from the settings stored in Cell 3
    n_monomers = get_param_from_state(wm, "parameters", "target_degree_polymerization")

    if not traj_path or not os.path.exists(traj_path):
        print("❌ Error: Trajectory file not found.")
        return None

    # 2. Load Trajectory
    print(f"   -> Loading trajectory: {os.path.basename(traj_path)}...")
    traj = md.load(traj_path, top=top_path)
    
    # ---------------------------------------------------------
    # CALCULATION 1: Radius of Gyration (Rg)
    # ---------------------------------------------------------
    rg = md.compute_rg(traj)
    avg_rg = np.mean(rg)
    std_rg = np.std(rg) # Standard Deviation (Fluctuation)
    
    # ---------------------------------------------------------
    # CALCULATION 2: End-to-End Distance (Ree)
    # ---------------------------------------------------------
    # Distances between first (0) and last atom (n_atoms-1)
    pairs = np.array([[0, traj.n_atoms - 1]])
    ree = md.compute_distances(traj, pairs)[:, 0] 
    avg_ree = np.mean(ree)
    std_ree = np.std(ree)

    # ---------------------------------------------------------
    # SAVE RESULTS TO JSON
    # ---------------------------------------------------------
    wm.save_result("radius_of_gyration", {
        "mean": round(float(avg_rg), 4),
        "std": round(float(std_rg), 4),
        "units": "nm"
    })
    
    wm.save_result("end_to_end_distance", {
        "mean": round(float(avg_ree), 4),
        "std": round(float(std_ree), 4),
        "units": "nm"
    })

    if n_monomers:
        wm.save_result("degree_of_polymerization", int(n_monomers))

    print(f"✅ Calculations Complete.")
    print(f"   Monomers: {n_monomers}")
    print(f"   Rg:  {avg_rg:.3f} ± {std_rg:.3f} nm")
    print(f"   Ree: {avg_ree:.3f} ± {std_ree:.3f} nm")

    # ---------------------------------------------------------
    # PLOTTING
    # ---------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    # --- Plot Rg ---
    ax[0].plot(traj.time, rg, color='#00796B', linewidth=1.5, label='$R_g$')
    
    # FIX: Added 'r' prefix to string for \sigma
    ax[0].fill_between(traj.time, avg_rg - std_rg, avg_rg + std_rg, color='#00796B', alpha=0.2, label=r'$\sigma$ (Fluctuation)')
    
    ax[0].axhline(avg_rg, color='black', linestyle='--', alpha=0.8, label=f'Mean: {avg_rg:.2f}')
    ax[0].set_title("Compactness ($R_g$)", fontsize=14)
    ax[0].set_xlabel("Time (ps)")
    ax[0].set_ylabel("Radius (nm)")
    ax[0].legend(loc='upper right')
    ax[0].grid(True, alpha=0.3)

    # --- Plot End-to-End ---
    ax[1].plot(traj.time, ree, color='#512DA8', linewidth=1.5, label='$R_{ee}$')
    
    # FIX: Added 'r' prefix to string for \sigma
    ax[1].fill_between(traj.time, avg_ree - std_ree, avg_ree + std_ree, color='#512DA8', alpha=0.2, label=r'$\sigma$ (Fluctuation)')
    
    ax[1].axhline(avg_ree, color='black', linestyle='--', alpha=0.8, label=f'Mean: {avg_ree:.2f}')
    ax[1].set_title("Extension ($R_{ee}$)", fontsize=14)
    ax[1].set_xlabel("Time (ps)")
    ax[1].set_ylabel("Distance (nm)")
    ax[1].legend(loc='upper right')
    ax[1].grid(True, alpha=0.3)

    plt.tight_layout()
    return fig
