import os
import sys
import time
import openmm as mm
from openmm import app
from openmm import unit

def run_vacuum_simulation(gro_file, top_file, output_prefix, temp_k, n_steps):
    """
    Runs a Vacuum NVT simulation using Gromacs input files.
    """
    
    # Define output filenames
    output_min = f"{output_prefix}_minimized.pdb"
    output_dcd = f"{output_prefix}_trajectory.dcd"
    output_final = f"{output_prefix}_final.pdb"
    
    print(f"   -> Loading Gromacs files...")
    print(f"      • GRO: {os.path.basename(gro_file)}")
    print(f"      • TOP: {os.path.basename(top_file)}")

    # 1. Load Gromacs Topology & Coordinates
    # We assume any included .itp files are in the same folder as the .top file
    include_dir = os.path.dirname(top_file)
    gro = app.GromacsGroFile(gro_file)
    top = app.GromacsTopFile(top_file, 
                             periodicBoxVectors=gro.getPeriodicBoxVectors(), 
                             includeDir=include_dir)

    # 2. Prepare for Vacuum
    # CRITICAL: We must remove the periodic box data from the topology.
    # If we don't, 'NoCutoff' will crash because it sees box vectors.
    top.topology.setUnitCellDimensions(None)

    # 3. Create OpenMM System
    print("   -> Creating System (Vacuum / NoCutoff)...")
    system = top.createSystem(nonbondedMethod=app.NoCutoff,
                              constraints=app.HBonds,
                              removeCMMotion=True)

    # 4. Integrator (Langevin NVT)
    dt_ps = 0.002  # 2 femtoseconds
    integrator = mm.LangevinMiddleIntegrator(temp_k * unit.kelvin,
                                             1.0 / unit.picosecond,
                                             dt_ps * unit.picoseconds)

    # 5. Select Hardware (GPU/CPU)
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
        props = {'Precision': 'mixed'}
    except:
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
            props = {'Precision': 'mixed'}
        except:
            platform = mm.Platform.getPlatformByName('CPU')
            props = {}
    
    print(f"   -> Using Platform: {platform.getName()}")

    # 6. Initialize Simulation
    simulation = app.Simulation(top.topology, system, integrator, platform, props)
    simulation.context.setPositions(gro.positions)

    # 7. Energy Minimization
    print("   -> Minimizing Energy...")
    simulation.minimizeEnergy()
    
    # Save Minimized PDB
    with open(output_min, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, 
                              simulation.context.getState(getPositions=True).getPositions(), 
                              f)
    print(f"      Saved: {output_min}")

    # 8. Run MD Simulation
    print(f"   -> Running MD ({n_steps} steps at {temp_k}K)...")
    
    # Add Reporters
    # DCD for trajectory (every 100 steps)
    simulation.reporters.append(app.DCDReporter(output_dcd, 100))
    # StateData for console logs (every 1000 steps)
    simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, 
                                                      potentialEnergy=True, temperature=True, 
                                                      speed=True, totalSteps=n_steps))

    start_time = time.time()
    simulation.step(n_steps)
    end_time = time.time()

    # 9. Save Final Topology/Frame
    # We save this so visualization tools can load the DCD topology correctly
    with open(output_final, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, 
                              simulation.context.getState(getPositions=True).getPositions(), 
                              f)
    print(f"      Saved: {output_final}")

    # Stats
    elapsed = end_time - start_time
    total_ns = (n_steps * dt_ps) / 1000.0
    ns_per_day = (total_ns / elapsed) * 86400 if elapsed > 0 else 0
    print(f"   -> Performance: {ns_per_day:.2f} ns/day")
    
    return output_dcd, output_final
