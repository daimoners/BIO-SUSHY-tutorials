"""
md_engine.py

OpenMM engines used by the BIO-SUSHY tutorial.

This module contains two complementary workflows:

1. run_vacuum_simulation()
   Short single-chain NVT simulation in vacuum.

2. run_bulk_annealing()
   Periodic NPT bulk annealing protocol.

The bulk protocol follows the logic:

high temperature / elevated pressure
→ cooling under pressure
→ final equilibration at room temperature and ambient pressure

The protocol is intentionally short for tutorial purposes. It should not be
interpreted as a production-grade polymer bulk equilibration.
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

    _save_pdb(simulation, output_min, enforce_periodic_box=False)

    sim_time_ps = int(n_steps) * dt_ps
    sim_time_ns = sim_time_ps / 1000.0

    print(f"   -> Running MD: {int(n_steps)} steps ({sim_time_ns:.3f} ns) at {temp_k} K...")

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
            totalSteps=int(n_steps),
        )
    )

    start_time = time.time()
    simulation.step(int(n_steps))
    end_time = time.time()

    _save_pdb(simulation, output_final, enforce_periodic_box=False)

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
#        TUNABLE BULK NPT PROTOCOL
# ==========================================

def build_bulk_protocol(
    hot_temperature_k=600.0,
    hot_pressure_bar=50.0,
    hot_steps=5000,
    cooling_start_temperature_k=600.0,
    cooling_end_temperature_k=300.0,
    cooling_pressure_bar=20.0,
    cooling_n_stages=4,
    cooling_steps_per_stage=3000,
    final_temperature_k=300.0,
    final_pressure_bar=1.0,
    final_steps=5000,
):
    """
    Build a tunable bulk NPT annealing protocol.

    Protocol:
    1. High-temperature / elevated-pressure annealing
    2. Stepwise cooling under pressure
    3. Final equilibration at room temperature and ambient pressure
    """
    cooling_n_stages = int(cooling_n_stages)
    if cooling_n_stages < 1:
        raise ValueError("cooling_n_stages must be at least 1.")

    stages = []

    stages.append({
        "name": "annealing_high_T_high_P",
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
        "name": "equilibration_300K_1bar",
        "temperature_k": float(final_temperature_k),
        "pressure_bar": float(final_pressure_bar),
        "steps": int(final_steps),
    })

    return stages


FAST_BULK_STAGES = build_bulk_protocol(
    hot_temperature_k=600.0,
    hot_pressure_bar=50.0,
    hot_steps=2000,
    cooling_start_temperature_k=600.0,
    cooling_end_temperature_k=300.0,
    cooling_pressure_bar=20.0,
    cooling_n_stages=2,
    cooling_steps_per_stage=1500,
    final_temperature_k=300.0,
    final_pressure_bar=1.0,
    final_steps=2000,
)


DEFAULT_BULK_STAGES = build_bulk_protocol(
    hot_temperature_k=600.0,
    hot_pressure_bar=50.0,
    hot_steps=5000,
    cooling_start_temperature_k=600.0,
    cooling_end_temperature_k=300.0,
    cooling_pressure_bar=20.0,
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
    """Select the best available OpenMM platform."""
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
    """Return the length of an OpenMM Vec3 or Quantity(Vec3) in nm."""
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
    """Return periodic box vector lengths in nm."""
    vectors = gro.getPeriodicBoxVectors()
    if vectors is None:
        raise ValueError("Bulk simulations require periodic box vectors in the GRO file.")

    return [_vec3_length_nm(v) for v in vectors]


def _safe_cutoff_nm(gro, requested_cutoff_nm):
    """Choose a nonbonded cutoff safely below half the smallest initial box length."""
    lengths = _box_lengths_nm(gro)
    min_box = min(lengths)

    max_allowed = 0.40 * min_box
    cutoff = min(float(requested_cutoff_nm), max_allowed)

    if cutoff < 0.25:
        raise ValueError(
            "The box is too small for a stable periodic nonbonded cutoff. "
            f"Minimum box length = {min_box:.3f} nm, suggested cutoff = {cutoff:.3f} nm. "
            "Build a larger initial bulk box or use fewer chains with the tutorial settings."
        )

    return float(cutoff)


def _normalise_stage(stage, index):
    """Accept both explicit and short stage key names."""
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


def _read_last_value_from_csv(csv_file, keyword):
    """Read the last value of a column containing keyword from a CSV file."""
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
            if keyword.lower() in key.lower():
                return float(value)
    except Exception:
        return None

    return None


def _save_pdb(simulation, filename, enforce_periodic_box=True):
    """Save current simulation positions to PDB."""
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)

    state = simulation.context.getState(
        getPositions=True,
        enforcePeriodicBox=bool(enforce_periodic_box),
    )

    with filename.open("w") as handle:
        app.PDBFile.writeFile(
            simulation.topology,
            state.getPositions(),
            handle,
        )

    return str(filename)


def _save_state_xml(simulation, filename):
    """Save current OpenMM state as XML."""
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
    """Print a human-readable summary of the NPT protocol."""
    print("\nBulk NPT protocol:")
    for i, stage in enumerate(stages, start=1):
        print(
            f"   {i:02d}. {stage['name']:<28s} "
            f"T = {stage['temperature_k']:7.1f} K | "
            f"P = {stage['pressure_bar']:7.1f} bar | "
            f"steps = {stage['steps']}"
        )


# ==========================================
#        BULK NPT ENGINE
# ==========================================

def run_bulk_annealing(
    gro_file,
    top_file,
    output_prefix,
    stages=None,
    timestep_ps=0.002,
    report_interval_steps=500,
    nonbonded_cutoff_nm=0.6,
    pressure_coupling_frequency=50,
    barostat_frequency=None,
    friction_per_ps=1.0,
    use_pme=True,
):
    """
    Run one continuous NPT bulk annealing/cooling simulation.

    The simulation uses one OpenMM Simulation object, one DCD trajectory and one
    CSV log. Temperature and pressure targets are updated stage-by-stage inside
    the same simulation.

    Notes
    -----
    Atomic motion is propagated by molecular dynamics.  The pressure-coupling
    mechanism changes the periodic box during the NPT simulation.
    """
    if stages is None:
        stages = DEFAULT_BULK_STAGES

    # Backward compatibility with older notebook cells.
    if barostat_frequency is not None:
        pressure_coupling_frequency = barostat_frequency

    stages = [_normalise_stage(stage, i) for i, stage in enumerate(stages)]

    timestep_ps = float(timestep_ps)
    if timestep_ps <= 0:
        raise ValueError("timestep_ps must be positive.")

    report_interval_steps = int(report_interval_steps)
    if report_interval_steps <= 0:
        raise ValueError("report_interval_steps must be positive.")

    pressure_coupling_frequency = int(pressure_coupling_frequency)
    if pressure_coupling_frequency <= 0:
        raise ValueError("pressure_coupling_frequency must be positive.")

    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    print("--- 🧱 Starting Bulk NPT MD Annealing ---")
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
    print(f"   -> Pressure coupling frequency: {pressure_coupling_frequency} steps")

    system = top.createSystem(
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
        removeCMMotion=True,
    )

    first_stage = stages[0]

    pressure_coupler = mm.MonteCarloBarostat(
        first_stage["pressure_bar"] * unit.bar,
        first_stage["temperature_k"] * unit.kelvin,
        pressure_coupling_frequency,
    )
    system.addForce(pressure_coupler)

    integrator = mm.LangevinMiddleIntegrator(
        first_stage["temperature_k"] * unit.kelvin,
        friction_per_ps / unit.picosecond,
        timestep_ps * unit.picoseconds,
    )

    platform, props = _choose_platform()
    print(f"   -> Using Platform: {platform.getName()}")

    simulation = app.Simulation(top.topology, system, integrator, platform, props)
    simulation.context.setPositions(gro.positions)
    simulation.context.setPeriodicBoxVectors(*gro.getPeriodicBoxVectors())

    print("   -> Minimizing Energy...")
    simulation.minimizeEnergy()

    minimized_pdb = _save_pdb(
        simulation,
        f"{output_prefix}_minimized.pdb",
        enforce_periodic_box=True,
    )

    dcd_file = f"{output_prefix}_trajectory.dcd"
    csv_file = f"{output_prefix}_log.csv"
    final_pdb_file = f"{output_prefix}_final.pdb"

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
        )
    )

    stage_summaries = []
    total_steps = 0
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

        integrator.setTemperature(temperature_k * unit.kelvin)

        # Update pressure-coupling targets.
        try:
            pressure_coupler.setDefaultTemperature(temperature_k * unit.kelvin)
            pressure_coupler.setDefaultPressure(pressure_bar * unit.bar)

            simulation.context.setParameter(
                mm.MonteCarloBarostat.Temperature(),
                temperature_k * unit.kelvin,
            )
            simulation.context.setParameter(
                mm.MonteCarloBarostat.Pressure(),
                pressure_bar * unit.bar,
            )
        except Exception:
            print("   ⚠️ Could not update pressure-coupling parameters dynamically.")

        simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin)

        stage_start_step = total_steps
        stage_start_wall = time.time()

        try:
            simulation.step(steps)
        except Exception as exc:
            raise RuntimeError(
                "Bulk NPT simulation failed during stage "
                f"'{stage_name}'. If the error mentions the nonbonded cutoff, "
                "reduce nonbonded_cutoff_nm, reduce pressure, or increase the "
                "initial box generated by bulk_builder.py."
            ) from exc

        stage_wall = time.time() - stage_start_wall
        total_steps += steps

        stage_summaries.append({
            "name": stage_name,
            "temperature_k": temperature_k,
            "pressure_bar": pressure_bar,
            "steps": steps,
            "start_step": int(stage_start_step),
            "end_step": int(total_steps),
            "wall_seconds": round(stage_wall, 3),
        })

    total_wall = time.time() - start_wall

    final_pdb = _save_pdb(
        simulation,
        final_pdb_file,
        enforce_periodic_box=True,
    )

    final_state_xml = _save_state_xml(simulation, f"{output_prefix}_final_state.xml")

    last_volume = _read_last_value_from_csv(csv_file, "Volume")
    last_density = _read_last_value_from_csv(csv_file, "Density")

    summary = {
        "status": "success",
        "minimized_pdb": os.path.abspath(minimized_pdb),
        "trajectory": os.path.abspath(dcd_file),
        "log_csv": os.path.abspath(csv_file),
        "final_pdb": os.path.abspath(final_pdb),
        "final_state_xml": os.path.abspath(final_state_xml),
        "stages": stage_summaries,
        "total_steps": int(total_steps),
        "timestep_ps": timestep_ps,
        "total_time_ps": float(total_steps * timestep_ps),
        "report_interval_steps": report_interval_steps,
        "nonbonded_cutoff_nm": cutoff_nm,
        "pressure_coupling_frequency": pressure_coupling_frequency,
        "friction_per_ps": float(friction_per_ps),
        "use_pme": bool(use_pme),
        "last_reported_volume": last_volume,
        "last_reported_density": last_density,
        "wall_seconds": round(total_wall, 3),
        "note": (
            "Tutorial NPT MD run. Density is logged for diagnostics but is not "
            "reported as a converged material property."
        ),
    }

    summary_file = f"{output_prefix}_summary.json"
    with open(summary_file, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=4)

    summary["summary_file"] = os.path.abspath(summary_file)

    print("\n" + "=" * 60)
    print("📊 MINI-BULK NPT MD REPORT")
    print(f"   Total simulated time: {summary['total_time_ps']:.2f} ps")
    print(f"   Trajectory:           {dcd_file}")
    print(f"   Log:                  {csv_file}")
    print(f"   Final structure:      {final_pdb}")
    print(f"   Wall time:            {total_wall:.2f} s")
    print("=" * 60 + "\n")

    return summary
