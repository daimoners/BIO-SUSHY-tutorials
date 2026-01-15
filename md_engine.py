import openmm as mm
from openmm import app
from openmm import unit
import sys
import os
import time  # <--- Added for timing

def run_simulation(pdb_file, xml_file, temp_k=300, press_bar=1, steps=50000, output_name="traj"):
    """
    Runs an NPT molecular dynamics simulation using OpenMM.
    Prints performance stats (Elapsed time, ns/day) at the end.
    """
    
    print(f"--- 🚀 Starting MD Simulation: {output_name} ---")
    
    # 1. Load Topology and Forcefield
    pdb = app.PDBFile(pdb_file)
    forcefield = app.ForceField(xml_file)
    
    # 2. Create System
    # NonbondedMethod=PME implies periodic boundary conditions are used
    system = forcefield.createSystem(pdb.topology, 
                                     nonbondedMethod=app.PME, 
                                     nonbondedCutoff=1.0*unit.nanometer, 
                                     constraints=app.HBonds)
    
    # 3. Integrator (Langevin - controls Temp)
    dt_ps = 0.002 # 2 fs timestep
    integrator = mm.LangevinMiddleIntegrator(temp_k*unit.kelvin, 
                                             1.0/unit.picosecond, 
                                             dt_ps*unit.picoseconds)
    
    # 4. Barostat (MonteCarlo - controls Pressure)
    system.addForce(mm.MonteCarloBarostat(press_bar*unit.bar, temp_k*unit.kelvin))
    
    # 5. Simulation Object
    # Try to use GPU (CUDA/OpenCL) if available, else CPU
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
            
    simulation = app.Simulation(pdb.topology, system, integrator, platform, props)
    simulation.context.setPositions(pdb.positions)
    
    # 6. Minimize Energy first (important to fix bad contacts)
    print("   -> Minimizing energy...")
    simulation.minimizeEnergy()
    
    # 7. Reporters (Output data)
    # Save trajectory every 1000 steps
    dcd_reporter = app.DCDReporter(f"{output_name}.dcd", 1000)
    
    # Print progress to console every 1000 steps
    data_reporter = app.StateDataReporter(sys.stdout, 1000, step=True, 
                                          potentialEnergy=True, temperature=True, 
                                          volume=True, speed=True)
    
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(data_reporter)
    
    # 8. Run Simulation (With Timer)
    print(f"   -> Running {steps} steps on {platform.getName()}...")
    
    start_time = time.time()  # <--- Start Timer
    simulation.step(steps)
    end_time = time.time()    # <--- Stop Timer
    
    # 9. Calculate Statistics
    elapsed_seconds = end_time - start_time
    total_time_ns = (steps * dt_ps) / 1000.0  # Convert steps*ps to ns
    ns_per_day = 0
    if elapsed_seconds > 0:
        ns_per_day = total_time_ns / (elapsed_seconds / 86400.0)
    
    print("-" * 40)
    print("📊 SIMULATION REPORT")
    print(f"   • Total Steps:      {steps}")
    print(f"   • Sim Time:         {total_time_ns:.3f} ns")
    print(f"   • Elapsed Time:     {elapsed_seconds:.2f} s ({elapsed_seconds/60:.2f} min)")
    print(f"   • Performance:      {ns_per_day:.2f} ns/day")
    print("-" * 40)
    
    # 10. Save Final State
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(f"{output_name}_final.pdb", 'w'))
    
    return f"{output_name}.dcd"
