%%writefile dft_runner.py
import os
import sys
import numpy as np
from IPython.display import clear_output

# Try to import PySCF, install if missing
try:
    import pyscf
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyscf"])
    import pyscf

from pyscf import gto, dft

def run_dft_validation(geometry_xyz, output_dir, basis="sto-3g"):
    """
    Runs FAST DFT (PBE/STO-3G with Density Fitting).
    Includes a live custom progress tracker.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = {"status": "failed"}
    
    # 1. Prepare Molecule
    atom_string = []
    with open(geometry_xyz, 'r') as f:
        lines = f.readlines()
        # Skip header (2 lines in XYZ)
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                atom_string.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}")
    mol_str = "; ".join(atom_string)

    try:
        # Build Molecule (Silent)
        mol = gto.M(atom=mol_str, basis=basis, verbose=0)
        mol.build()

        # 2. Setup FAST DFT
        # RKS = Restricted Kohn-Sham (Standard DFT)
        # .density_fit() = Approximates integrals (10x Speedup)
        # xc='pbe' = Pure functional (Faster than Hybrid B3LYP)
        mf = dft.RKS(mol).density_fit()
        mf.xc = 'pbe' 
        
        # 3. Define Live Progress Tracker
        print("⏳ Initializing PySCF (DFT)...")
        
        # This function runs inside the DFT loop every iteration
        def progress_tracker(envs):
            cycle = envs.get('cycle', 0)
            e_tot = envs.get('e_tot', 0)
            
            # Update screen IN-PLACE
            clear_output(wait=True)
            print(f"⚛️  DFT Simulation Running (Fast Mode)...\n   ► Step:   {cycle+1}\n   ► Energy: {e_tot:.5f} Eh")
            
        mf.callback = progress_tracker

        # 4. Run Simulation
        mf.kernel()
        
        # Clear the "Running" message when done
        clear_output(wait=True)

        # 5. Extract Properties
        if mf.converged:
            # Dipole
            dip_vec = mf.dip_moment(unit='Debye', verbose=0)
            dipole = np.linalg.norm(dip_vec)

            # Gap (Approximate)
            if mf.mo_occ is not None:
                lumo_idx = np.where(mf.mo_occ == 0)[0][0]
                homo_idx = lumo_idx - 1
                homo_ev = mf.mo_energy[homo_idx] * 27.2114
                lumo_ev = mf.mo_energy[lumo_idx] * 27.2114
                gap = abs(lumo_ev - homo_ev)
            else:
                gap = 0.0

            results = {
                "status": "success",
                "dipole": dipole,
                "band_gap": gap,
                "energy": mf.e_tot,
                "converged": True
            }
        else:
            results["error"] = "DFT SCF did not converge."

    except Exception as e:
        results["error"] = str(e)
        
    return results
