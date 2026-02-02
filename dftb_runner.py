import os
import json
import numpy as np
from ase import io, units
from ase.md.langevin import Langevin
from ase.calculators.calculator import Calculator, all_changes
from xtb.ase.calculator import XTB # Requires xtb-python

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    """
    Runs a semi-empirical Quantum MD simulation using GFN2-xTB.
    
    Args:
        input_pdb (str): Path to the starting structure (from Classical MD).
        output_prefix (str): Base name for outputs.
        temp_k (float): Temperature in Kelvin.
        steps (int): Number of MD steps (typically 0.5 - 1.0 fs per step).
    
    Returns:
        dict: Electronic properties (HOMO-LUMO gap, Dipole).
    """
    print(f"--- ⚛️ Starting GFN2-xTB (DFTB) Simulation ---")
    print(f"   • Input: {os.path.basename(input_pdb)}")
    
    # 1. Load Structure into ASE
    atoms = io.read(input_pdb)
    
    # 2. Attach XTB Calculator (GFN2-xTB is the default method)
    # accuracy=1.0 is standard. Lower is faster but less accurate.
    atoms.calc = XTB(method="GFN2-xTB", accuracy=1.0, temperature=temp_k)
    
    # 3. Setup MD (Langevin Dynamics)
    # Time step: 1.0 fs (Standard for XTB involving H atoms)
    timestep = 1.0 * units.fs
    friction = 0.02  # Friction coefficient
    
    dyn = Langevin(atoms, timestep=timestep, temperature_K=temp_k, friction=friction,
                   trajectory=f"{output_prefix}_traj.traj",
                   logfile=f"{output_prefix}.log")
    
    # 4. Run MD
    print(f"   -> Running {steps} steps of Quantum MD...")
    dyn.run(steps)
    
    # 5. Final Analysis (Single Point on final frame)
    # We force a calculation to get electronic properties
    energy = atoms.get_potential_energy()
    gap = atoms.calc.get_homo_lumo_gap()
    dipole = np.linalg.norm(atoms.get_dipole_moment())
    
    # Save Final PDB
    final_pdb = f"{output_prefix}_final.pdb"
    io.write(final_pdb, atoms)
    
    results = {
        "final_energy_eV": energy,
        "homo_lumo_gap_eV": gap,
        "dipole_moment_debye": dipole
    }
    
    print(f"   ✅ Done. Gap: {gap:.3f} eV | Dipole: {dipole:.3f} D")
    return results, final_pdb
