# BIO-SUSHY-tutorial
CNR physics-based modelling team - Modelling Tutorial Hub

**A modular, Python-based Computational Tutorial for simulating and analyzing polymer dynamics.**

![Status](https://img.shields.io/badge/Status-Active-success)
![Platform](https://img.shields.io/badge/Platform-Google%20Colab%20%7C%20Linux-orange)
![License](https://img.shields.io/badge/License-MIT-blue)

## 📖 Overview

This project provides an automated workflow to generate, simulate, and analyze synthetic polymers. Designed to be modular, it allows users to input a monomer (SMILES), build a polymer chain, run Molecular Dynamics (MD) simulations in a vacuum, and visualize the trajectory interactively.

It features a robust **State Management System** that tracks files and results automatically, ensuring reproducible research.

## ✨ Key Features

* **🧬 Automatic Structure Generation:** Converts monomer SMILES strings into 3D PDB structures (both homo- and co-polymers) and scales them to a specified Degree of Polymerization (DP).
* **⚡ OpenMM Simulations:** Runs vacuum MD simulations using the Langevin integrator.
* **📊 Automated Analysis:** Calculates Radius of Gyration ($R_g$) and End-to-End Distance ($R_{end}$) over time.
* **🎥 Interactive Visualization:** Includes a 3D viewer (py3Dmol) with animation controls, a "step counter" HUD, and a VMD-style frame inspector.
* **💾 State Persistence:** Uses a JSON-based workflow manager (`state.json`) to save paths, parameters, and results automatically.

---

## 📂 Project Structure

The project is designed with modular Python scripts to separate concerns:

| File | Description |
| :--- | :--- |
| **`workflow.py`** | **The Brain.** Manages `state.json` and `results.json`, ensuring data flows correctly between modules. |
| **`polymer_builder.py`** | Handles geometry generation. Converts SMILES to PDB and replicates units to form a polymer chain. |
| **`ff_builder.py`** | **Force Field Application.** Applies OPLS-AA parameters to the generated structure using Foyer. |
| **`md_engine.py`** | **The MD Engine.** Loads the system, sets up the OpenMM integrator (Langevin), and runs the vacuum simulation. |
| **`analysis.py`** | **Post-processing.** Calculates physical metrics ($R_g$, $R_{end}$) and generates plots with error analysis. |
| **`visualiser.py`** | **3D Graphics.** Contains logic for the interactive widget, style settings, and the animation player. |
| **`oplsaa.xml`** | Contains the specific force field parameters used by `ff_builder.py`. |

---

## 🔬 How It Works

### 1. The Workflow Manager (`workflow.py`)
Instead of passing variables manually between notebook cells, the `WorkflowManager` saves the current state of the experiment to `state.json`.
* *Example:* The **Builder** saves the PDB path to the state. The **Engine** reads that path from the state automatically.

### 2. Simulation Details
* **Force Field:** OPLS-AA (via `oplsaa.xml` and Foyer).
* **Integrator:** Langevin Integrator (Simulating thermal bath collisions).
* **Conditions:** Vacuum simulation (No explicit solvent, Infinite Cutoff).

### 3. Visualization
The viewer (`visualiser.py`) uses **py3Dmol** but is enhanced with:
* **Custom Styling:** Standard CPK coloring with optimized sphere/stick ratios.
* **HUD:** A Heads-Up Display showing the current simulation step.
* **Frame Inspector:** A widget to scrub through individual frames for detailed inspection.

---

## 📊 Outputs

The simulation generates a folder named after your polymer (e.g., `PolyStyrene/`). Inside, you will find:

* `MD/polymer_vac_trajectory.dcd`: The raw motion data of the atoms.
* `MD/polymer_vac_final.pdb`: The final equilibrated structure.
* `results.json`: A summary of the calculated metrics (e.g., Average $R_g$, Simulation Time).
* `workflow.log`: A text log of the operations performed.

---

## 🤝 Contributing

Contributions are welcome! Please fork the repository and submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
