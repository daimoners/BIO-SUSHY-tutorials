import os
import sys
import numpy as np

# Try to import PySCF, install if missing
try:
    import pyscf
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyscf"])
    import pyscf

from pyscf import gto, dft, lib

def run_dft_validation(geometry_xyz, output_dir, basis="sto-3g"):
    """
    Runs DFT with streaming output to the notebook.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = {"status": "failed"}
    
    # 1. Stream PySCF output to screen
    # We create a logger that prints to sys.stdout
    print(f"   ► Running PySCF (DFT) in: {output_dir}")
    print("   ► Iteration logs will appear below:")
    print("-" * 50)

    try:
        # 2. Parse XYZ
        atom_string = []
        with open(geometry_xyz, 'r') as f:
            lines = f.readlines()
            for line in lines[2:]:
                parts = line.split()
                if len(parts) >= 4:
                    atom_string.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}")
        mol_str = "; ".join(atom_string)

        # 3. Build Molecule
        mol = gto.M(atom=mol_str, basis=basis, verbose=4) # Verbose 4 = Print SCF iterations
        mol.output = None # None means print to stdout (screen)
        mol.build()

        # 4. Run DFT
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        # This kernel call will now print the "cycle 1 E=..." lines automatically
        mf.kernel()

        print("-" * 50) # End of log

        # 5. Extract Properties
        if mf.converged:
            # Dipole
            dip_vec = mf.dip_moment(unit='Debye', verbose=0)
            dipole = np.linalg.norm(dip_vec)

            # Gap
            lumo_idx = np.where(mf.mo_occ == 0)[0][0]
            homo_idx = lumo_idx - 1
            homo_ev = mf.mo_energy[homo_idx] * 27.2114
            lumo_ev = mf.mo_energy[lumo_idx] * 27.2114
            gap = abs(lumo_ev - homo_ev)

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
