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

from pyscf import gto, dft, data

def run_dft_validation(geometry_xyz, output_dir, basis="sto-3g", xc="pbe", properties=None):
    """
    Runs DFT with User-Selected Parameters.
    Args:
        basis (str): Basis set (e.g., 'sto-3g', '6-31g').
        xc (str): Exchange-Correlation functional (e.g., 'pbe', 'b3lyp').
        properties (list): List of metrics to compute ['dipole', 'gap', 'energy'].
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
        # Build Molecule
        # verbose=0 suppresses internal printing (we use our own streamer)
        mol = gto.M(atom=mol_str, basis=basis, verbose=0)
        mol.build()

        # 2. Setup DFT (Customizable)
        # Density Fitting is ALWAYS used here to keep the tutorial fast
        mf = dft.RKS(mol).density_fit()
        mf.xc = xc 
        
        # 3. Live Progress Tracker
        print(f"⏳ Initializing PySCF ({xc}/{basis})...")
        
        def progress_tracker(envs):
            cycle = envs.get('cycle', 0)
            e_tot = envs.get('e_tot', 0)
            
            # Format energy
            e_str = f"{e_tot:.4f}" if abs(e_tot) > 0.1 else "..."
            
            clear_output(wait=True)
            print(f"⚛️  DFT Simulation Running...\n   ► Setup:  {xc.upper()} / {basis}\n   ► Step:   {cycle+1}\n   ► Energy: {e_str} Eh")
            
        mf.callback = progress_tracker

        # 4. Run Kernel
        mf.kernel()
        clear_output(wait=True)

        if mf.converged:
            # --- CALCULATE METRICS ---
            calc_data = {"converged": True, "energy": mf.e_tot}
            
            # 1. Dipole Moment
            if "dipole" in properties:
                dip_vec = mf.dip_moment(unit='Debye', verbose=0)
                calc_data["dipole"] = np.linalg.norm(dip_vec)
                
            # 2. HOMO-LUMO Gap
            if "gap" in properties:
                if mf.mo_occ is not None:
                    # Find HOMO (highest occupied) and LUMO (lowest unoccupied)
                    lumo_idx = np.where(mf.mo_occ == 0)[0][0]
                    homo_idx = lumo_idx - 1
                    
                    homo_ev = mf.mo_energy[homo_idx] * 27.2114
                    lumo_ev = mf.mo_energy[lumo_idx] * 27.2114
                    
                    calc_data["gap"] = abs(lumo_ev - homo_ev)
                    calc_data["homo"] = homo_ev
                    calc_data["lumo"] = lumo_ev
                else:
                    calc_data["gap"] = 0.0

            # 3. Mulliken Charges (Optional advanced metric)
            if "charges" in properties:
                # Simple population analysis
                (pop, chg) = mf.mulliken_pop(verbose=0)
                calc_data["max_charge"] = np.max(np.abs(chg)) # Max local charge

            results = {**calc_data, "status": "success"}
        else:
            results["error"] = "DFT SCF did not converge."

    except Exception as e:
        results["error"] = str(e)
        
    return results
