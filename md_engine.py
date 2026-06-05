"""
md_engine.py

OpenMM engines used by the BIO-SUSHY tutorial.

This module contains two complementary workflows:

1. run_vacuum_simulation()
   Short single-chain NVT simulation in vacuum.

2. run_bulk_annealing()
   Periodic NPT mini-bulk annealing protocol.

The mini-bulk protocol follows the logic:

high temperature / high pressure
→ cooling under pressure
→ final equilibration at room temperature and ambient pressure

The protocol is intentionally short for tutorial purposes.
It should not be interpreted as a production-grade polymer bulk equilibration.
"""

import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app
from openmm import unit


# ==========================================
#        SINGLE-CHAIN VACUUM MD
# ==========================================

def run_vacuum_simulation(
    gro_file,
    top_file,
    output_prefix,
    temp_k,
    n_steps,
    timestep_ps=0.002,
    report_interval_steps=100,
):
    """
    Runs a vacuum NVT simulation using Gromacs input files.

    Parameters
    ----------
    gro_file : str
        Input GRO coordinate file.
    top_file : str
        Input Gromacs topology file.
    output_prefix : str
        Prefix for output files.
    temp_k : float
        Simulation temperature in K.
    n_steps : int
        Number of MD steps.
    timestep_ps : float
        MD timestep in ps. Default 0.002 ps = 2 fs.
    report_interval_steps : int
        Number of steps between saved trajectory frames.

    Returns
    -------
    tuple
        output_dcd, output_final_pdb
    """

    output_min = f"{output_prefix}_minimized.pdb"
    output_dcd = f"{output_prefix}_trajectory.dcd"
    output_final = f"{output_prefix}_final.pdb"

    print("--- 🌪️ Starting Vacuum Simulation ---")
    print(f"   • GRO: {os.path.basename(gro_file)}")
    print(f"   • TOP: {os.path.basename(top_file)}")

    include_dir = os.path.dirname(top_file)

    gro = app.GromacsGroFile(gro_file)
    top = app.GromacsTopFile(
        top_file,
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=include_dir,
    )

    # Vacuum simulation: remove periodicity and use NoCutoff.
    top.topology.setUnitCellDimensions(None)

    print("   -> Creating System (Vacuum / NoCutoff)...")
    system = top.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        removeCMMotion=True,
    )

    dt_ps = float(timestep_ps)
    if dt_ps <= 0:
        raise ValueError("timestep_ps must be positive.")

    integrator = mm.LangevinMiddleIntegrator(
        temp_k * unit.kelvin,
        1.0 / unit.picosecond,
        dt_ps * unit.picoseconds,
    )

    platform, props = _choose_platform()
    print(f"   -> Using Platform: {platform.getName()}")

    simulation = app.Simulation(top.topology, system, integrator, platform, props)
    simulation.context.setPositions(gro.positions)

    print("   -> Minimizing Energy...")
    simulation.minimizeEnergy()

    with open(output_min, "w") as handle:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True).getPositions(),
            handle,
        )

    sim_time_ps = n_steps * dt_ps
    sim_time_ns = sim_time_ps / 1000.0

    print(f"   -> Running MD: {n_steps} steps ({sim_time_ns:.3f} ns) at {temp_k} K...")

    report_interval_steps = int(report_interval_steps)
    if report_interval_steps <= 0:
        raise ValueError("report_interval_steps must be positive.")

    print(f"   -> Saving trajectory frames every {report_interval_steps} steps.")

    simulation.reporters.append(app.DCDReporter(output_dcd, report_interval_steps))
    simulation.reporters.append(
        app.StateDataReporter(
            sys.stdout,
            1000,
            step=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
            remainingTime=True,
            totalSteps=n_steps,
        )
    )

    start_time = time.time()
    simulation.step(int(n_steps))
    end_time = time.time()

    with open(output_final, "w") as handle:
        app.PDBFile.writeFile(
            simulation.topology,
            simulation.context.getState(getPositions=True).getPositions(),
            handle,
        )

    elapsed_seconds = end_time - start_time
    ns_per_day = 0.0
    if elapsed_seconds > 0:
        ns_per_day = (sim_time_ns / elapsed_seconds) * 86400.0

    print("\n" + "=" * 40)
    print("📊 SIMULATION REPORT")
    print(f"   • Simulated Time:   {sim_time_ps:.1f} ps")
    print(f"   • Report Interval:  {report_interval_steps} steps")
    print(f"   • Wall Clock Time:  {elapsed_seconds:.2f} s ({elapsed_seconds / 60:.2f} min)")
    print(f"   • Performance:      {ns_per_day:.2f} ns/day")
    print("=" * 40 + "\n")

    return output_dcd, output_final


# ==========================================
#        TUNABLE MINI-BULK NPT PROTOCOL
# ==========================================

def build_bulk_protocol(
    hot_temperature_k=600.0,
    hot_pressure_bar=300.0,
    hot_steps=5000,
    cooling_start_temperature_k=600.0,
    cooling_end_temperature_k=300.0,
    cooling_pressure_bar=100.0,
    cooling_n_stages=4,
    cooling_steps_per_stage=3000,
    final_temperature_k=300.0,
    final_pressure_bar=1.0,
    final_steps=5000,
):
    """
    Build a tunable mini-bulk NPT annealing protocol.

    Protocol:
    1. High-temperature / high-pressure compression
    2. Stepwise cooling under high pressure
    3. Final equilibration at room temperature and ambient pressure

    Parameters
    ----------
    hot_temperature_k : float
        Temperature for the initial hot-compression stage.
    hot_pressure_bar : float
        Pressure for the initial hot-compression stage.
    hot_steps : int
        Number of steps for the hot-compression stage.
    cooling_start_temperature_k : float
        First temperature of the cooling ramp.
    cooling_end_temperature_k : float
        Last temperature of the cooling ramp.
    cooling_pressure_bar : float
        Pressure used during the cooling stages.
    cooling_n_stages : int
        Number of cooling stages.
    cooling_steps_per_stage : int
        Number of MD steps for each cooling stage.
    final_temperature_k : float
        Temperature for final equilibration.
    final_pressure_bar : float
        Pressure for final equilibration.
        For ambient pressure, use 1.0 bar.
    final_steps : int
        Number of MD steps for final equilibration.

    Returns
    -------
    list of dict
        Stage list compatible with run_bulk_annealing().
    """

    cooling_n_stages = int(cooling_n_stages)
    if cooling_n_stages < 1:
        raise ValueError("cooling_n_stages must be at least 1.")

    stages = []

    stages.append({
        "name": "hot_compression",
        "temperature_k": float(hot_temperature_k),
        "pressure_bar": float(hot_pressure_bar),
        "steps": int(hot_steps),
    })

    cooling_temperatures = np.linspace(
        float(cooling_start_temperature_k),
        float(cooling_end_temperature_k),
        cooling_n_stages,
    )

    for temp_k in cooling_temperatures:
        stages.append({
            "name": f"cooling_{int(round(temp_k))}K",
            "temperature_k": float(temp_k),
            "pressure_bar": float(cooling_pressure_bar),
            "steps": int(cooling_steps_per_stage),
        })

    stages.append({
        "name": "final_equilibration",
        "temperature_k": float(final_temperature_k),
        "pressure_bar": float(final_pressure_bar),
        "steps": int(final_steps),
    })

    return stages


FAST_BULK_STAGES = build_bulk_protocol(
    hot_temperature_k=600.0,
    hot_pressure_bar=300.0,
    hot_steps=2000,
    cooling_start_temperature_k=600.0,
    cooling_end_temperature_k=300.0,
    cooling_pressure_bar=100.0,
    cooling_n_stages=2,
    cooling_steps_per_stage=1500,
    final_temperature_k=300.0,
    final_pressure_bar=1.0,
    final_steps=2000,
)


DEFAULT_BULK_STAGES = build_bulk_protocol(
    hot_temperature_k=650.0,
    hot_pressure_bar=300.0,
    hot_steps=5000,
    cooling_start_temperature_k=600.0,
    cooling_end_temperature_k=300.0,
    cooling_pressure_bar=100.0,
    cooling_n_stages=4,
    cooling_steps_per_stage=3000,
    final_temperature_k=300.0,
    final_pressure_bar=1.0,
    final_steps=5000,
)


# ==========================================
#        INTERNAL HELPERS
# ==========================================

def _choose_platform():
    """
    Select the best available OpenMM platform.
    """
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


def _vec3_length_nm(vector):
    """
    Return the length of an OpenMM Vec3 or Quantity(Vec3) in nm.

    This avoids using vector.norm(), which is not always available for
    Quantity-wrapped Vec3 objects.
    """
    try:
        vec_nm = vector.value_in_unit(unit.nanometer)
        return float(np.sqrt(vec_nm[0] ** 2 + vec_nm[1] ** 2 + vec_nm[2] ** 2))
    except Exception:
        values = []
        for component in vector:
            try:
                values.append(component.value_in_unit(unit.nanometer))
            except Exception:
                values.append(float(component))
        return float(np.sqrt(sum(v * v for v in values)))


def _box_lengths_nm(gro):
    """
    Return box vector lengths in nm.
    """
    vectors = gro.getPeriodicBoxVectors()
    if vectors is None:
        raise ValueError("Bulk simulations require periodic box vectors in the GRO file.")

    return [_vec3_length_nm(v) for v in vectors]


def _safe_cutoff_nm(gro, requested_cutoff_nm):
    """
    Ensure the nonbonded cutoff is smaller than half the smallest box length.
    """
    lengths = _box_lengths_nm(gro)
    min_box = min(lengths)

    max_allowed = 0.45 * min_box
    cutoff = min(float(requested_cutoff_nm), max_allowed)

    if cutoff < 0.4:
        raise ValueError(
            "The box is too small for a stable periodic nonbonded cutoff. "
            f"Minimum box length = {min_box:.3f} nm, suggested cutoff = {cutoff:.3f} nm. "
            "Build a larger initial bulk box."
        )

    return cutoff


def _normalise_stage(stage, index):
    """
    Accept both explicit and short stage key names.
    """
    name = stage.get("name", f"stage_{index + 1}")

    temperature = stage.get("temperature_k", stage.get("temperature", None))
    pressure = stage.get("pressure_bar", stage.get("pressure", None))
    steps = stage.get("steps", None)

    if temperature is None or pressure is None or steps is None:
        raise ValueError(
            f"Invalid stage {index + 1}. "
            "Required keys: name, temperature_k, pressure_bar, steps."
        )

    return {
        "name": str(name),
        "temperature_k": float(temperature),
        "pressure_bar": float(pressure),
        "steps": int(steps),
    }


def _read_last_density_from_csv(csv_file):
    """
    Read the last density value from an OpenMM StateDataReporter CSV file.
    """
    path = Path(csv_file)
    if not path.exists():
        return None

    try:
        with path.open() as handle:
            rows = list(csv.DictReader(handle))

        if not rows:
            return None

        last_row = rows[-1]

        for key, value in last_row.items():
            if "Density" in key:
                return float(value)

    except Exception:
        return None

    return None


def _save_pdb(simulation, filename):
    """
    Save current simulation positions to PDB.
    """
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)

    state = simulation.context.getState(
        getPositions=True,
        enforcePeriodicBox=True,
    )

    with filename.open("w") as handle:
        app.PDBFile.writeFile(
            simulation.topology,
            state.getPositions(),
            handle,
        )

    return str(filename)


def _save_state_xml(simulation, filename):
    """
    Save current OpenMM state as XML.
    """
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)

    state = simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )

    with filename.open("w") as handle:
        handle.write(mm.XmlSerializer.serialize(state))

    return str(filename)


def _print_protocol(stages):
    """
    Print a human-readable summary of the NPT protocol.
    """
    print("\nMini-bulk NPT protocol:")
    for i, stage in enumerate(stages, start=1):
        print(
            f"   {i:02d}. {stage['name']:<24s} "
            f"T = {stage['temperature_k']:7.1f} K | "
            f"P = {stage['pressure_bar']:7.1f} bar | "
            f"steps = {stage['steps']}"
        )


# ==========================================
#        MINI-BULK NPT ENGINE
# ==========================================

def run_bulk_annealing(
    gro_file,
    top_file,
    output_prefix,
    stages=None,
    timestep_ps=0.002,
    report_interval_steps=500,
    nonbonded_cutoff_nm=1.0,
    barostat_frequency=25,
    friction_per_ps=1.0,
    use_pme=True,
):
    """
    Run a compact NPT mini-bulk annealing protocol.

    Parameters
    ----------
    gro_file : str
        Periodic GRO coordinate file for the mini-bulk.
    top_file : str
        Gromacs topology file for the mini-bulk.
    output_prefix : str
        Prefix for output files, e.g. "Polymer/Bulk/anneal".
    stages : list of dict
        List of dictionaries with keys:
        name, temperature_k, pressure_bar, steps.
        If None, DEFAULT_BULK_STAGES is used.
    timestep_ps : float
        Integration timestep in ps. Default 0.002 ps = 2 fs.
    report_interval_steps : int
        Reporting frequency for DCD and CSV files.
    nonbonded_cutoff_nm : float
        Requested nonbonded cutoff.
        It is automatically reduced if the box is small.
    barostat_frequency : int
        Frequency in MD steps for Monte Carlo volume moves.
    friction_per_ps : float
        Langevin friction coefficient in 1/ps.
    use_pme : bool
        If True, use PME. If False, use CutoffPeriodic.

    Returns
    -------
    dict
        Summary with output paths and protocol metadata.
    """

    if stages is None:
        stages = DEFAULT_BULK_STAGES

    stages = [_normalise_stage(stage, i) for i, stage in enumerate(stages)]

    timestep_ps = float(timestep_ps)
    if timestep_ps <= 0:
        raise ValueError("timestep_ps must be positive.")

    report_interval_steps = int(report_interval_steps)
    if report_interval_steps <= 0:
        raise ValueError("report_interval_steps must be positive.")

    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("--- 🧱 Starting Mini-Bulk NPT Annealing ---")
    print(f"   • GRO: {os.path.basename(gro_file)}")
    print(f"   • TOP: {os.path.basename(top_file)}")
    print(f"   • Stages: {len(stages)}")

    _print_protocol(stages)

    include_dir = os.path.dirname(top_file)

    gro = app.GromacsGroFile(gro_file)
    top = app.GromacsTopFile(
        top_file,
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=include_dir,
    )

    cutoff_nm = _safe_cutoff_nm(gro, nonbonded_cutoff_nm)
    nonbonded_method = app.PME if use_pme else app.CutoffPeriodic

    print("\n   -> Creating periodic OpenMM system...")
    print(f"   -> Nonbonded method: {'PME' if use_pme else 'CutoffPeriodic'}")
    print(f"   -> Nonbonded cutoff: {cutoff_nm:.3f} nm")

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
        timestep_ps * unit.picoseconds,
    )

    platform, props = _choose_platform()
    print(f"   -> Using Platform: {platform.getName()}")

    simulation = app.Simulation(
        top.topology,
        system,
        integrator,
        platform,
        props,
    )

    simulation.context.setPositions(gro.positions)
    simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())

    print("   -> Minimizing Energy...")
    simulation.minimizeEnergy()

    minimized_pdb = _save_pdb(
        simulation,
        f"{output_prefix}_minimized.pdb",
    )

    stage_logs = []
    stage_trajectories = []
    stage_final_pdbs = []
    stage_summaries = []

    start_wall = time.time()

    for i, stage in enumerate(stages, start=1):
        stage_name = stage["name"].replace(" ", "_")
        temperature_k = stage["temperature_k"]
        pressure_bar = stage["pressure_bar"]
        steps = int(stage["steps"])

        print("\n" + "=" * 60)
        print(f"Stage {i}/{len(stages)}: {stage_name}")
        print(f"   T = {temperature_k:.1f} K")
        print(f"   P = {pressure_bar:.1f} bar")
        print(f"   steps = {steps}")

        # Update thermostat temperature.
        integrator.setTemperature(temperature_k * unit.kelvin)

        # Update barostat target temperature and pressure.
        try:
            barostat.setDefaultTemperature(temperature_k * unit.kelvin)
            barostat.setDefaultPressure(pressure_bar * unit.bar)

            simulation.context.setParameter(
                mm.MonteCarloBarostat.Temperature(),
                temperature_k * unit.kelvin,
            )
            simulation.context.setParameter(
                mm.MonteCarloBarostat.Pressure(),
                pressure_bar * unit.bar,
            )
        except Exception:
            print(
                "   ⚠️ Could not update barostat parameters dynamically. "
                "Continuing with the existing barostat settings."
            )

        # Reassign velocities at each stage.
        simulation.context.setVelocitiesToTemperature(
            temperature_k * unit.kelvin,
        )

        dcd_file = f"{output_prefix}_{stage_name}.dcd"
        csv_file = f"{output_prefix}_{stage_name}.csv"
        final_pdb = f"{output_prefix}_{stage_name}_final.pdb"

        simulation.reporters = []

        simulation.reporters.append(
            app.DCDReporter(
                dcd_file,
                report_interval_steps,
            )
        )

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
        last_density = _read_last_density_from_csv(csv_file)

        stage_logs.append(csv_file)
        stage_trajectories.append(dcd_file)
        stage_final_pdbs.append(final_pdb)

        stage_summaries.append(
            {
                "name": stage_name,
                "temperature_k": temperature_k,
                "pressure_bar": pressure_bar,
                "steps": steps,
                "csv_log": os.path.abspath(csv_file),
                "trajectory": os.path.abspath(dcd_file),
                "final_pdb": os.path.abspath(final_pdb),
                "wall_seconds": round(stage_wall, 3),
                "last_reported_density": last_density,
            }
        )

    total_wall = time.time() - start_wall

    final_pdb = stage_final_pdbs[-1]
    final_state_xml = _save_state_xml(
        simulation,
        f"{output_prefix}_final_state.xml",
    )

    summary = {
        "status": "success",
        "minimized_pdb": os.path.abspath(minimized_pdb),
        "final_pdb": os.path.abspath(final_pdb),
        "final_state_xml": os.path.abspath(final_state_xml),
        "stage_logs": [os.path.abspath(x) for x in stage_logs],
        "stage_trajectories": [os.path.abspath(x) for x in stage_trajectories],
        "stage_final_pdbs": [os.path.abspath(x) for x in stage_final_pdbs],
        "stages": stage_summaries,
        "timestep_ps": timestep_ps,
        "report_interval_steps": report_interval_steps,
        "nonbonded_cutoff_nm": cutoff_nm,
        "barostat_frequency": int(barostat_frequency),
        "friction_per_ps": float(friction_per_ps),
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
