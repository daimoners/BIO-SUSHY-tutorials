import os
import sys
import numpy as np

# Try to import PySCF, install if missing (automating dependency)
try:
    import pyscf
except ImportError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyscf"])
    import pyscf

from pyscf import gto, dft

def run_dft_validation(geometry_xyz, output_dir, basis="sto-3g"):
    """
    Runs a Single Point DFT (B3LYP) calculation on the provided geometry.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = {"status": "failed"}
    
    try:
        # 1. Parse XYZ
        atom_string = []
        with open(geometry_xyz, 'r') as f:
            lines = f.readlines()
            # Skip header (atom count + comment)
            for line in lines[2:]:
                parts = line.split()
                if len(parts) >= 4:
                    atom_string.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}")
        
        mol_str = "; ".join(atom_string)

        # 2. Build Molecule
        mol = gto.M(atom=mol_str, basis=basis, verbose=0)
        mol.build()

        # 3. Run DFT
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        mf.kernel()

        # 4. Extract Properties
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
            "converged": mf.converged
        }

    except Exception as e:
        results["error"] = str(e)
        
    return results
