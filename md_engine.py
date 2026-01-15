# @title 8a. Install MD Engine Module 🛠️
# @markdown This creates the `md_engine.py` script which handles Minimization and NVT dynamics.

import os

engine_code = 
import os
import sys
import openmm
import openmm.app as app
import openmm.unit as unit
from parmed import load_file

def run_vacuum_simulation(gro_file, top_file, output_prefix="sim", 
                          temp_k=300, n_steps=10000, report_interval=100):
    """
    Runs a Vacuum NVT simulation:
    1. Loads OPLS topology
    2. Minimizes Energy (Steepest Descent)
    3. Runs Dynamics (Langevin)
    """
    print(f"🔧 MD ENGINE: Initializing {output_prefix}...")
    
    # 1. Load Data
    if not os.path.exists(gro_file) or not os.path.exists(top_file):
        raise FileNotFoundError(f"Missing input files: {gro_file} or {top_file}")

    print("   -> Loading structure (ParmEd)...")
    structure = load_file(top_file, xyz=gro_file)

    # 2. Build System (Vacuum = No Cutoff)
    system = structure.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        implicitSolvent=None
    )

    # 3. Integrator (Thermostat)
    integrator = openmm.LangevinIntegrator(
        temp_k * unit.kelvin, 
        1.0 / unit.picosecond,  # Friction
        2.0 * unit.femtoseconds # Step size
    )

    # 4. Simulation Context
    platform = openmm.Platform.getPlatformByName('CPU') # Safe default
    simulation = app.Simulation(structure.topology, system, integrator, platform)
    simulation.context.setPositions(structure.positions)

    # 5. Energy Minimization
    print(f"   -> 📉 Minimizing Energy (Relaxing structure)...")
    print(f"      Initial Energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}")
    simulation.minimizeEnergy()
    min_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"      Final Energy:   {min_energy}")
    
    # Save Minimized State
    min_pdb = f"{output_prefix}_minimized.pdb"
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(min_pdb, 'w'))
    print(f"      Saved minimized structure: {min_pdb}")

    # 6. Production Run
    print(f"   -> 🚀 Starting NVT Simulation ({n_steps} steps at {temp_k}K)...")
    
    # Output Files
    traj_file = f"{output_prefix}_trajectory.dcd"
    
    # Reporters
    simulation.reporters.append(app.DCDReporter(traj_file, report_interval))
    simulation.reporters.append(app.StateDataReporter(
        sys.stdout, report_interval*10, step=True, potentialEnergy=True, temperature=True
    ))

    # Run
    simulation.step(n_steps)
    
    # Save Final PDB for visualization reference
    final_pdb = f"{output_prefix}_final.pdb"
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(final_pdb, 'w'))
    
    print(f"✅ Simulation Complete. Trajectory: {traj_file}")
    return traj_file, final_pdb
"""

# Write the script to the BIO-SUSHY folder so Python finds it
repo_path = os.path.abspath("BIO-SUSHY-tutorials")
if not os.path.exists(repo_path): os.makedirs(repo_path)

with open(os.path.join(repo_path, "md_engine.py"), "w") as f:
    f.write(engine_code)

print("✅ MD Engine Installed.")
