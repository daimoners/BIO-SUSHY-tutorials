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

def run_dft_validation(geometry_xyz, output_dir, basis="sto-3g", xc="pbe", properties=None):
    """
    Runs DFT Single Point Calculation.
    """
    if properties is None: properties = ["dipole", "gap"]
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = {"status": "failed"}
    
    # 1. Prepare Molecule
    atom_string = []
    with open(geometry_xyz, 'r') as f:
        lines = f.readlines()
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 4:
                atom_string.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}")
    mol_str = "; ".join(atom_string)

    try:
        # Build Molecule (Silent)
        mol = gto.M(atom=mol_str, basis=basis, verbose=0)
        mol.build()

        # 2. Setup DFT
        mf = dft.RKS(mol).density_fit()
        mf.xc = xc 
        
        # 3. Live Progress Tracker
        print(f"⏳ Initializing PySCF ({xc.upper()}/{basis})...")
        
        def progress_tracker(envs):
            cycle = envs.get('cycle', 0)
            e_tot = envs.get('e_tot', 0)
            # Update screen IN-PLACE
            clear_output(wait=True)
            print(f"⚛️  DFT Single Point Calculation...\n   ► Geometry: Fixed (Pre-optimized)\n   ► SCF Cycle: {cycle+1}\n   ► Energy:   {e_tot:.5f} Eh")
            
        mf.callback = progress_tracker

        # 4. Run Kernel
        mf.kernel()
        clear_output(wait=True)

        # 5. Extract Properties
        if mf.converged:
            calc_data = {"converged": True, "energy": mf.e_tot}
            
            if "dipole" in properties:
                dip_vec = mf.dip_moment(unit='Debye', verbose=0)
                calc_data["dipole"] = np.linalg.norm(dip_vec)
            
            if "gap" in properties:
                if mf.mo_occ is not None:
                    lumo_idx = np.where(mf.mo_occ == 0)[0][0]
                    homo_idx = lumo_idx - 1
                    homo_ev = mf.mo_energy[homo_idx] * 27.2114
                    lumo_ev = mf.mo_energy[lumo_idx] * 27.2114
                    calc_data["gap"] = abs(lumo_ev - homo_ev)
                else:
                    calc_data["gap"] = 0.0
            
            results = {**calc_data, "status": "success"}
        else:
            results["error"] = "DFT SCF did not converge."

    except Exception as e:
        results["error"] = str(e)
        
    return results
