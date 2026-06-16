# BIO-SUSHY-tutorial

CNR DAIMON team - Modelling Tutorial Hub

**A modular, Python-based computational tutorial for building, simulating, visualising and analysing polymer systems.**

![Status](https://img.shields.io/badge/Status-Active-success)
![Platform](https://img.shields.io/badge/Platform-Google%20Colab%20%7C%20Linux-orange)
![License](https://img.shields.io/badge/License-GPLv3-blue)

---

## 📖 Overview

`BIO-SUSHY-tutorial` provides a Colab-ready workflow for early-stage, physics-based polymer modelling. Starting from monomer or comonomer SMILES strings, the workflow can:

1. build and relax a single polymer chain;
2. assign OPLS-AA force-field parameters;
3. run a short single-chain MD simulation in vacuum;
4. extract conformational descriptors such as radius of gyration and end-to-end distance;
5. build a demonstrative small bulk system from multiple polymer chains;
6. run a short periodic NPT annealing / cooling protocol;
7. analyse qualitative bulk relaxation descriptors such as box volume and mass density;
8. visualise static structures and trajectories, including periodic boxes when available.

NOTE: The tutorial is designed for workflow demonstration, training and descriptor generation. It is **not** a validated production protocol for converged polymer properties, toxicity, persistence, exposure, LCA or full SSbD assessment.

---

## ✨ Key Features

- **SMILES-to-polymer construction**: builds homopolymers and copolymers from user-provided repeat-unit SMILES.
- **Monomer inspection and structure relaxation**: generates 3D monomer views and minimised polymer-chain PDB structures.
- **OPLS-AA parameterisation**: applies OPLS-AA atom typing through Foyer and writes Gromacs-compatible `.gro` and `.top` files.
- **Single-chain vacuum MD**: runs short OpenMM NVT simulations for conformational sampling.
- **bulk construction**: replicates one parameterised chain into a cubic periodic box with automatically estimated starting box size.
- **Bulk NPT annealing**: supports tunable high-temperature / high-pressure annealing, cooling and final equilibration stages.
- **Automated analysis**: saves single-chain and bulk descriptors to `results.json`.
- **Interactive visualisation**: uses `py3Dmol` for static structures, trajectories and optional periodic-box rendering.
- **Workflow state management**: uses `state.json`, `results.json` and `workflow.log` to keep notebook cells modular and reproducible.

---

## 🔁 Workflow at a Glance

```text
Input SMILES, polymer type, degree of polymerisation, stoichiometry
        │
        ▼
polymer_builder.py
Build / inspect / minimise single polymer chain
        │
        ▼
ff_builder.py
Apply OPLS-AA parameters and write GRO/TOP files
        │
        ├──────────────► Single-chain route
        │                 md_engine.run_vacuum_simulation()
        │                 polymer_analyzer.analyze_trajectory()
        │                 Rg, Ree, analysis window
        │
        └──────────────► Bulk route
                          bulk_builder.build_polymer_bulk()
                          md_engine.run_bulk_annealing()
                          bulk_analyzer.analyze_bulk_density()
                          volume, density, relaxation plots
```

---

## 📂 Project Structure

| File | Role |
| :--- | :--- |
| `workflow.py` | Workflow state manager. Creates the project folder and manages `state.json`, `results.json` and `workflow.log`. |
| `polymer_builder.py` | Builds polymer chains from SMILES strings, supports homopolymers and copolymers, generates 3D structures and performs molecular-mechanics relaxation. |
| `ff_builder.py` | Applies OPLS-AA parameters using Foyer and writes Gromacs-compatible `conf.gro` and `topol.top` files. |
| `bulk_builder.py` | Builds a tutorial bulk system from one parameterised polymer chain. The user chooses the number of chains; the initial box is estimated automatically. |
| `md_engine.py` | Contains OpenMM engines for both single-chain vacuum NVT simulations and periodic bulk NPT annealing. |
| `polymer_analyzer.py` | Analyses single-chain trajectories, computing radius of gyration, end-to-end distance and trajectory-window metadata. |
| `bulk_analyzer.py` | Analyses bulk OpenMM CSV logs, plotting box-volume and density relaxation and saving final-window descriptors. |
| `visualiser.py` | Provides `py3Dmol` static and trajectory visualisation utilities, including GRO-to-PDB conversion and periodic-box drawing. |
| `oplsaa.xml` | OPLS-AA force-field XML used by `ff_builder.py`. |
| `LICENSE` | Project license. |

---

## 🧩 Main Modules

### `workflow.py`

The `WorkflowManager` keeps the notebook modular by storing paths, parameters and results on disk. Typical files written inside each polymer project directory are:

- `state.json`: paths and workflow parameters;
- `results.json`: computed descriptors and metadata;
- `workflow.log`: chronological log of workflow operations.

Typical methods include:

```python
wm.add_path("gro_file", gro_path)
wm.get_path("trajectory")
wm.update_state("parameters", {...})
wm.save_result("radius_of_gyration", {...})
```

### `polymer_builder.py`

Handles molecular construction before force-field assignment.

Main user-facing functions:

```python
inspect_monomer(smiles, label="monomer")
build_polymer(smiles_A, polymer_type, n_monomers, smiles_B="", ratio=0.5, name="polymer")
minimize_polymer(input_pdb, name="polymer")
get_relaxed_files()
```

Supported polymer types include:

- `homopolymer`
- `copolymer-alternate`
- `copolymer-block`
- `copolymer-random`

### `ff_builder.py`

Applies OPLS-AA atom types and writes the parameterised simulation files.

Main user-facing functions:

```python
build_opls_system(pdb_file, output_name="system_opls")
inspect_system(gro_file, top_file)
```

The module also sanitises the local `oplsaa.xml` file for OpenMM compatibility and includes helper functions for charge handling.

### `md_engine.py`

Contains two complementary OpenMM workflows.

#### Single-chain vacuum simulation

```python
run_vacuum_simulation(
    gro_file,
    top_file,
    output_prefix,
    temp_k,
    n_steps,
    timestep_ps=0.002,
    report_interval_steps=100,
)
```

This runs a short NVT simulation in vacuum using a Langevin integrator and `NoCutoff` nonbonded treatment. It returns the DCD trajectory and final PDB structure.

#### bulk NPT annealing

```python
build_bulk_protocol(...)
run_bulk_annealing(
    gro_file,
    top_file,
    output_prefix,
    stages=None,
    timestep_ps=0.002,
    report_interval_steps=500,
    nonbonded_cutoff_nm=0.6,
    pressure_coupling_frequency=50,
    friction_per_ps=1.0,
    use_pme=True,
)
```

The default bulk protocol follows:

```text
high-temperature / elevated-pressure annealing
→ stepwise cooling under pressure
→ final equilibration at room temperature and ambient pressure
```

Two pre-defined stage lists are available:

- `FAST_BULK_STAGES`: shorter demonstration protocol;
- `DEFAULT_BULK_STAGES`: longer tutorial protocol.

### `bulk_builder.py`

Builds a periodic bulk starting configuration from a single parameterised polymer chain.

```python
build_polymer_bulk(
    single_chain_gro,
    single_chain_top,
    n_chains=6,
    output_dir="Bulk",
    random_seed=42,
    initial_packing_density_kg_m3=180.0,
    minimum_box_nm=4.5,
    padding_nm=0.8,
)
```

Returns:

```python
bulk_gro, bulk_top, metadata
```

The initial box size is calculated internally from chain mass, chain extent, number of chains and a dilute packing-density estimate. This initial density is only used to create a stable starting geometry and should not be interpreted as a predicted material density.

### `polymer_analyzer.py`

Analyses single-chain trajectories.

```python
analyze_trajectory(
    wm,
    equilibration_fraction=0.5,
    use_heavy_atom_endpoints=True,
)
```

Computed descriptors include:

- radius of gyration, `radius_of_gyration`;
- end-to-end distance, `end_to_end_distance`;
- analysis-window metadata;
- degree of polymerisation, when available in the workflow state.

By default, statistics are computed on the second half of the saved trajectory to avoid mixing the initial relaxation period with the analysis window.

### `bulk_analyzer.py`

Analyses OpenMM CSV logs from the bulk NPT simulation.

```python
analyze_bulk_density(
    wm=None,
    stage_logs=None,
    final_fraction=0.5,
    save_results=True,
)

print_bulk_interpretation()
```

Computed descriptors include:

- final-window mass density, `bulk_density`, in `g/cm^3`;
- final-window box volume, `bulk_volume`, in `nm^3`;
- density and volume relaxation plots.

These values are tutorial descriptors, not converged production-grade bulk properties.

### `visualiser.py`

Provides unified interactive visualisation utilities.

```python
show_structure(file_path_or_string, width=800, height=400, show_box="auto")
show_trajectory(pdb_content, width=800, height=400, show_box=False)
```

Supported static structure inputs include:

- `.mol`
- `.pdb`
- `.gro`
- raw PDB / MOL strings

For `.gro` files, the module can convert the structure to a minimal PDB block and draw the orthorhombic simulation box.

---

## 📊 Outputs

A typical workflow creates a project folder named after the polymer, for example:

```text
PolymerName/
├── state.json
├── results.json
├── workflow.log
├── *_polymer.pdb
├── *_relaxed.pdb
├── system_opls/
│   ├── conf.gro
│   └── topol.top
├── MD/
│   ├── polymer_vac_minimized.pdb
│   ├── polymer_vac_trajectory.dcd
│   └── polymer_vac_final.pdb
└── Bulk/
    ├── bulk.gro
    ├── topol.top
    ├── bulk_*_minimized.pdb
    ├── bulk_*_trajectory.dcd
    ├── bulk_*_log.csv
    ├── bulk_*_final.pdb
    ├── bulk_*_final_state.xml
    └── bulk_*_summary.json
```

Exact names may differ depending on the notebook `output_prefix` values.

---

## 🌱 Descriptor Scope and SSbD Context

The descriptors generated here can support early-stage, SSbD-inspired discussion by linking polymer structure, simulation protocol and qualitative material-level behaviour. Examples include:

- single-chain compactness through `Rg`;
- single-chain extension through `Ree`;
- qualitative bulk compaction through volume relaxation;
- qualitative bulk density trends from short NPT runs.

However, this tutorial does **not** provide a complete Safe-and-Sustainable-by-Design assessment. It does not evaluate hazard, toxicity, release, persistence, exposure, end-of-life, life-cycle impacts or regulatory compliance.

---

## ⚠️ Tutorial Limitations

- Simulations are intentionally short for training and Colab compatibility.
- bulk systems are small and should be treated as demonstrative.
- Density values from short NPT annealing are qualitative diagnostics, not validated material predictions.
- Force-field assignment is automated and should be inspected for unusual chemistries.
- Results should not be compared directly with experiment unless the protocol has been validated for the specific polymer class.

---

## 📄 License

This project is licensed under the GNU General Public License v3.0. See the `LICENSE` file for details.
