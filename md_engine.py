"""
md_engine.py

OpenMM engines used by the BIO-SUSHY tutorial.

This module contains two complementary workflows:

1. run_vacuum_simulation
   A short single-chain NVT simulation in vacuum, used as the beginner
   demonstration.

2. run_bulk_annealing
   A compact periodic NPT mini-bulk annealing protocol, used to demonstrate
   how the same workflow can move from an isolated chain to material-like
   descriptors such as density and volume relaxation.

The mini-bulk protocol is intentionally short and pedagogical. It should not
be interpreted as a production-grade equilibration protocol for quantitative
polymer properties.
"""

import csv
import json
import os
import sys
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import openmm as mm
from openmm import app
from openmm import unit

def run_vacuum_simulation(gro_file, top_file, output_prefix, temp_k, n_steps, timestep_ps=0.002, report_interval_steps=100):
    """
    Runs a Vacuum NVT simulation using Gromacs input files.
    Includes performance benchmarking and reporting.

    Parameters
    ----------
    timestep_ps : float, optional
        MD timestep in picoseconds. The default is 0.002 ps, i.e. 2 fs.
    report_interval_steps : int, optional
        Number of MD integration steps between saved trajectory frames.
        The default is 100, matching the original tutorial behaviour.
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
    dt_ps = float(timestep_ps)  # ps; default 0.002 ps = 2 femtoseconds
    if dt_ps <= 0:
        raise ValueError("timestep_ps must be a positive number.")
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
    sim_time_ns = sim_time_ps / 1000.0

    print(f"   -> Running MD: {n_steps} steps ({sim_time_ns:.3f} ns) at {temp_k}K...")
    
    # Reporters
    report_interval_steps = int(report_interval_steps)
    if report_interval_steps <= 0:
        raise ValueError("report_interval_steps must be a positive integer.")
    print(f"   -> Saving trajectory frames every {report_interval_steps} steps.")
    simulation.reporters.append(app.DCDReporter(output_dcd, report_interval_steps))
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
    print(f"   • Simulated Time:   {int(sim_time_ps)} ps")
    print(f"   • Report Interval:  {report_interval_steps} steps")
    print(f"   • Wall Clock Time:  {elapsed_seconds:.2f} s ({elapsed_seconds/60:.2f} min)")
    print(f"   • Performance:      {ns_per_day:.2f} ns/day")
    print("="*40 + "\n")
    
    return output_dcd, output_final


# ==========================================
#        MINI-BULK NPT ENGINE
# ==========================================

DEFAULT_BULK_STAGES = [
    {"name": "hot_compression", "temperature_k": 600.0, "pressure_bar": 300.0, "steps": 5000},
    {"name": "cool_500K", "temperature_k": 500.0, "pressure_bar": 100.0, "steps": 3000},
    {"name": "cool_400K", "temperature_k": 400.0, "pressure_bar": 50.0, "steps": 3000},
    {"name": "final_300K", "temperature_k": 300.0, "pressure_bar": 1.0, "steps": 5000},
]

FAST_BULK_STAGES = [
    {"name": "hot_compression", "temperature_k": 600.0, "pressure_bar": 300.0, "steps": 2000},
    {"name": "cool_450K", "temperature_k": 450.0, "pressure_bar": 100.0, "steps": 1500},
    {"name": "final_300K", "temperature_k": 300.0, "pressure_bar": 1.0, "steps": 2000},
]


def _choose_platform():
    """Return the best available OpenMM platform and default properties."""
    try:
        platform = mm.Platform.getPlatformByName("CUDA")
        props = {"Precision": "mixed"}
    except Exception:
        try:
            platform = mm.Platform.getPlatformByName("OpenCL")
            props = {"Precision": "mixed"}
        except Exception:
            platform = mm.Platform.getPlatformByName("CPU")
            props = {}
    return platform, props


def _box_min_length_nm(gro: app.GromacsGroFile) -> float:
    vectors = gro.getPeriodicBoxVectors()
    if vectors is None:
        raise ValueError("Bulk simulations require periodic box vectors in the GRO file.")
    lengths = []
    for v in vectors:
        lengths.append(v.norm().value_in_unit(unit.nanometer))
    return float(min(lengths))


def _safe_cutoff_nm(gro: app.GromacsGroFile, requested_cutoff_nm: float) -> float:
    """Ensure the nonbonded cutoff is smaller than half the smallest box length."""
    min_box = _box_min_length_nm(gro)
    max_allowed = 0.45 * min_box
    cutoff = min(float(requested_cutoff_nm), max_allowed)
    if cutoff < 0.4:
        raise ValueError(
            f"The box is too small for a stable periodic nonbonded cutoff: "
            f"min box length={min_box:.3f} nm, suggested cutoff={cutoff:.3f} nm. "
            "Build a larger initial bulk box."
        )
    return cutoff


def _normalise_stage(stage: dict, index: int) -> dict:
    """Accept both short and explicit stage key names."""
    name = stage.get("name", f"stage_{index+1}")
    temperature = stage.get("temperature_k", stage.get("temperature", None))
    pressure = stage.get("pressure_bar", stage.get("pressure", None))
    steps = stage.get("steps", None)
    if temperature is None or pressure is None or steps is None:
        raise ValueError(
            f"Invalid stage {index+1}. Required keys: name, temperature_k, pressure_bar, steps."
        )
    return {
        "name": str(name),
        "temperature_k": float(temperature),
        "pressure_bar": float(pressure),
        "steps": int(steps),
    }


def _read_last_density_from_csv(csv_file: str | os.PathLike) -> float | None:
    """Read the last density value from an OpenMM StateDataReporter CSV, if present."""
    path = Path(csv_file)
    if not path.exists():
        return None
    try:
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            return None
        last = rows[-1]
        for key in last:
            if "Density" in key:
                return float(last[key])
    except Exception:
        return None
    return None


def _save_pdb(simulation: app.Simulation, filename: str | os.PathLike) -> str:
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with filename.open("w") as handle:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle)
    return str(filename)


def run_bulk_annealing(
    gro_file: str,
    top_file: str,
    output_prefix: str,
    stages: Sequence[dict] | None = None,
    timestep_ps: float = 0.002,
    report_interval_steps: int = 500,
    nonbonded_cutoff_nm: float = 1.0,
    barostat_frequency: int = 25,
    friction_per_ps: float = 1.0,
    use_pme: bool = True,
) -> dict:
    """
    Run a compact NPT mini-bulk annealing protocol.

    Parameters
    ----------
    gro_file, top_file
        Periodic Gromacs coordinate/topology files for the mini-bulk.
    output_prefix
        Prefix for output files, e.g. "MyPolymer/Bulk/anneal".
    stages
        List of dictionaries with keys: name, temperature_k, pressure_bar, steps.
        If None, DEFAULT_BULK_STAGES is used.
    timestep_ps
        Integration timestep in ps. Default 0.002 ps = 2 fs.
    report_interval_steps
        Reporting frequency for DCD and CSV files.
    nonbonded_cutoff_nm
        Requested nonbonded cutoff. It will be reduced automatically if the box
        is too small.
    barostat_frequency
        Frequency in MD steps for Monte Carlo volume moves.
    friction_per_ps
        Langevin friction coefficient in 1/ps.
    use_pme
        If True, use PME. If False, use CutoffPeriodic.

    Returns
    -------
    dict
        Paths and metadata for downstream analysis.
    """
    stages = [_normalise_stage(s, i) for i, s in enumerate(stages or DEFAULT_BULK_STAGES)]
    report_interval_steps = int(report_interval_steps)
    if report_interval_steps <= 0:
        raise ValueError("report_interval_steps must be positive.")
    if timestep_ps <= 0:
        raise ValueError("timestep_ps must be positive.")

    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("--- 🧱 Starting Mini-Bulk NPT Annealing ---")
    print(f"   • GRO: {os.path.basename(gro_file)}")
    print(f"   • TOP: {os.path.basename(top_file)}")
    print(f"   • Stages: {len(stages)}")

    include_dir = os.path.dirname(top_file)
    gro = app.GromacsGroFile(gro_file)
    top = app.GromacsTopFile(
        top_file,
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=include_dir,
    )

    cutoff_nm = _safe_cutoff_nm(gro, nonbonded_cutoff_nm)
    nonbonded_method = app.PME if use_pme else app.CutoffPeriodic

    print("   -> Creating periodic OpenMM system...")
    system = top.createSystem(
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
        removeCMMotion=True,
    )

    first_stage = stages[0]
    barostat = mm.MonteCarloBarostat(
        first_stage["pressure_bar"] * unit.bar,
        first_stage["temperature_k"] * unit.kelvin,
        int(barostat_frequency),
    )
    system.addForce(barostat)

    integrator = mm.LangevinMiddleIntegrator(
        first_stage["temperature_k"] * unit.kelvin,
        friction_per_ps / unit.picosecond,
        float(timestep_ps) * unit.picoseconds,
    )

    platform, props = _choose_platform()
    print(f"   -> Using Platform: {platform.getName()}")
    print(f"   -> Nonbonded: {'PME' if use_pme else 'CutoffPeriodic'}; cutoff={cutoff_nm:.3f} nm")

    simulation = app.Simulation(top.topology, system, integrator, platform, props)
    simulation.context.setPositions(gro.positions)
    simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())

    print("   -> Minimizing Energy...")
    simulation.minimizeEnergy()
    minimized_pdb = _save_pdb(simulation, f"{output_prefix}_minimized.pdb")

    stage_logs = []
    stage_dcds = []
    stage_final_pdbs = []
    stage_summaries = []
    start_wall = time.time()

    for i, stage in enumerate(stages, start=1):
        name = stage["name"].replace(" ", "_")
        temperature_k = stage["temperature_k"]
        pressure_bar = stage["pressure_bar"]
        steps = stage["steps"]

        print("\n" + "=" * 60)
        print(f"Stage {i}/{len(stages)}: {name}")
        print(f"   T = {temperature_k:.1f} K | P = {pressure_bar:.1f} bar | steps = {steps}")

        integrator.setTemperature(temperature_k * unit.kelvin)
        try:
            simulation.context.setParameter(mm.MonteCarloBarostat.Temperature(), temperature_k * unit.kelvin)
            simulation.context.setParameter(mm.MonteCarloBarostat.Pressure(), pressure_bar * unit.bar)
        except Exception:
            # Older OpenMM versions may be stricter. The default values still
            # make the first stage valid; later stages remain demonstrative.
            print("   ⚠️ Could not update barostat parameters dynamically on this OpenMM version.")

        simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)

        dcd_file = f"{output_prefix}_{name}.dcd"
        csv_file = f"{output_prefix}_{name}.csv"
        final_pdb = f"{output_prefix}_{name}_final.pdb"

        simulation.reporters = []
        simulation.reporters.append(app.DCDReporter(dcd_file, report_interval_steps))
        simulation.reporters.append(
            app.StateDataReporter(
                csv_file,
                report_interval_steps,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
                separator=",",
            )
        )
        simulation.reporters.append(
            app.StateDataReporter(
                sys.stdout,
                max(report_interval_steps, 1000),
                step=True,
                temperature=True,
                volume=True,
                density=True,
                speed=True,
                totalSteps=steps,
                remainingTime=True,
            )
        )

        stage_start = time.time()
        simulation.step(steps)
        stage_wall = time.time() - stage_start

        final_pdb = _save_pdb(simulation, final_pdb)
        density = _read_last_density_from_csv(csv_file)

        stage_logs.append(csv_file)
        stage_dcds.append(dcd_file)
        stage_final_pdbs.append(final_pdb)
        stage_summaries.append(
            {
                **stage,
                "csv_log": os.path.abspath(csv_file),
                "trajectory": os.path.abspath(dcd_file),
                "final_pdb": os.path.abspath(final_pdb),
                "wall_seconds": round(stage_wall, 3),
                "last_reported_density": density,
            }
        )

    total_wall = time.time() - start_wall
    final_pdb = stage_final_pdbs[-1]
    final_state_xml = f"{output_prefix}_final_state.xml"
    with open(final_state_xml, "w") as handle:
        handle.write(mm.XmlSerializer.serialize(simulation.context.getState(getPositions=True, getVelocities=True)))

    summary = {
        "status": "success",
        "minimized_pdb": os.path.abspath(minimized_pdb),
        "final_pdb": os.path.abspath(final_pdb),
        "final_state_xml": os.path.abspath(final_state_xml),
        "stage_logs": [os.path.abspath(x) for x in stage_logs],
        "stage_trajectories": [os.path.abspath(x) for x in stage_dcds],
        "stage_final_pdbs": [os.path.abspath(x) for x in stage_final_pdbs],
        "stages": stage_summaries,
        "timestep_ps": float(timestep_ps),
        "report_interval_steps": int(report_interval_steps),
        "nonbonded_cutoff_nm": float(cutoff_nm),
        "use_pme": bool(use_pme),
        "wall_seconds": round(total_wall, 3),
    }

    summary_file = f"{output_prefix}_summary.json"
    with open(summary_file, "w") as handle:
        json.dump(summary, handle, indent=4)
    summary["summary_file"] = os.path.abspath(summary_file)

    print("\n" + "=" * 60)
    print("📊 MINI-BULK ANNEALING REPORT")
    print(f"   Final structure: {final_pdb}")
    print(f"   Logs:            {len(stage_logs)} CSV files")
    print(f"   Wall time:       {total_wall:.2f} s")
    print("=" * 60 + "\n")

    return summary
