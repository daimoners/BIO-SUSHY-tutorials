import os
import sys
import time
import openmm as mm
from openmm import app
from openmm import unit

def run_vacuum_simulation(gro_file, top_file, output_prefix, temp_k, n_steps):
    """
    Runs a Vacuum NVT simulation using Gromacs input files.
    Includes performance benchmarking and reporting.
    """
    
    # Define output filenames
    output_min = f"{output_prefix}_minimized.pdb"
    output_dcd = f"{output_prefix}_trajectory.dcd"
    output_final = f"{output_prefix}_final.pdb"
    
    print(f"--- 🌪️ Starting Vacuum Simulation ---")
    print(f"   • GRO: {os.path.basename(gro_file)}")
    print(f"   • TOP: {os.path.basename(top_file)}")

    # 1. Load Gromacs Topology & Coordinates
    # We assume any included .itp files are in the same folder as the .top file
    include_dir = os.path.dirname(top_file)
    gro = app.GromacsGroFile(gro_file)
    top = app.GromacsTopFile(top_file, 
                             periodicBoxVectors=gro.getPeriodicBoxVectors(), 
                             includeDir=include_dir)

    # 2. Prepare for Vacuum
    # CRITICAL: Remove periodic box to allow NoCutoff
    top.topology.setUnitCellDimensions(None)

    # 3. Create OpenMM System
    print("   -> Creating System (Vacuum / NoCutoff)...")
    system = top.createSystem(nonbondedMethod=app.NoCutoff,
                              constraints=app.HBonds,
                              removeCMMotion=True)

    # 4. Integrator (Langevin NVT)
    dt_ps = 0.002  # 2 femtoseconds timestep
    integrator = mm.LangevinMiddleIntegrator(temp_k * unit.kelvin,
                                             1.0 / unit.picosecond,
                                             dt_ps * unit.picoseconds)

    # 5. Select Hardware
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

    # 8. Run MD Simulation
    sim_time_ps = (n_steps * dt_ps) 
    sim_time_ns = tim_time_ps / 1000.0

    print(f"   -> Running MD: {n_steps} steps ({sim_time_ns:.3f} ns) at {temp_k}K...")
    
    # Reporters
    simulation.reporters.append(app.DCDReporter(output_dcd, 100))
    simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, 
                                                      potentialEnergy=True, temperature=True, 
                                                      speed=True, remainingTime=True, 
                                                      totalSteps=n_steps))

    # --- TIMER START ---
    start_time = time.time()
    
    simulation.step(n_steps)
    
    # --- TIMER END ---
    end_time = time.time()

    # 9. Save Final Topology/Frame
    with open(output_final, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, 
                              simulation.context.getState(getPositions=True).getPositions(), 
                              f)

    # 10. Performance Statistics
    elapsed_seconds = end_time - start_time
    ns_per_day = 0.0
    if elapsed_seconds > 0:
        ns_per_day = (sim_time_ns / elapsed_seconds) * 86400.0

    print("\n" + "="*40)
    print("📊 SIMULATION REPORT")
    print(f"   • Simulated Time:   {sim_time_ps:.4f} ps")
    print(f"   • Wall Clock Time:  {elapsed_seconds:.2f} s ({elapsed_seconds/60:.2f} min)")
    print(f"   • Performance:      {ns_per_day:.2f} ns/day")
    print("="*40 + "\n")
    
    return output_dcd, output_final
