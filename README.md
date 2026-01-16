# BIO-SUSHY-tutorial
CNR physics-based modelling team - Modelling Tutorial Hub

# 🧪 Polymer Dynamics Lab

**A modular, Python-based computational laboratory for simulating and analyzing polymer dynamics.**

![Status](https://img.shields.io/badge/Status-Active-success)
![Platform](https://img.shields.io/badge/Platform-Google%20Colab%20%7C%20Linux-orange)
![License](https://img.shields.io/badge/License-MIT-blue)

## 📖 Overview

This project provides a streamlined workflow to generate, simulate, and analyze synthetic polymers. Designed primarily for **Google Colab**, it allows users to input a monomer (SMILES), build a polymer chain, run Molecular Dynamics (MD) simulations in a vacuum, and visualize the trajectory interactively.

It features a robust **State Management System** that tracks files and results automatically, ensuring reproducible research.

## ✨ Key Features

* **🧬 Automatic Structure Generation:** Converts monomer SMILES strings into 3D PDB structures and scales them to a specified Degree of Polymerization (DP).
* **⚡ OpenMM Simulations:** Runs vacuum MD simulations using the Langevin integrator.
* **📊 Automated Analysis:** Calculates Radius of Gyration ($R_g$) and End-to-End Distance ($R_{end}$) over time.
* **🎥 Interactive Visualization:** Includes a 3D viewer (py3Dmol) with animation controls, a "step counter" HUD, and a VMD-style frame inspector.
* **💾 State Persistence:** Uses a JSON-based workflow manager (`state.json`) to save paths, parameters, and results automatically.

---

## 🚀 Quick Start (Google Colab)

The easiest way to run this lab is via Google Colab.

1.  **Clone this repository** into your Colab environment.
2.  **Open the Notebook** (e.g., `main.ipynb`).
3.  **Run the cells in order:**
    * The notebook will automatically install dependencies (`openmm`, `rdkit`, `mdtraj`, etc.) via `condacolab`.
    * Enter your Monomer SMILES (e.g., Styrene: `C=Cc1ccccc1`) and Chain Length.
    * Watch the simulation run and analyze the results.

---

## 📦 Installation (Local)

If you prefer to run this locally, you need **Conda** installed.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git](https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git)
    cd YOUR_REPO_NAME
    ```

2.  **Create the environment:**
    ```bash
    conda create -n polymerlab python=3.10
    conda activate polymerlab
    ```

3.  **Install dependencies:**
    ```bash
    conda install -c conda-forge openmm rdkit mdtraj py3dmol ipywidgets matplotlib
    ```

---

## 📂 Project Structure

The project is designed with modular Python scripts to separate concerns:

| File | Description |
| :--- | :--- |
| `builder.py` | Handles geometry generation. Converts SMILES to PDB and replicates units to form a polymer chain. |
| `simulator.py` | The MD Engine. Loads the PDB, sets up the OpenMM system (Forcefield, Integrator), and runs the simulation. |
| `analyzer.py` | Post-processing. Calculates physical metrics ($R_g$, $R_{end}$) from the trajectory files. |
| `visualiser.py` | 3D Graphics. Contains logic for the interactive widget, style settings, and the animation player. |
| `workflow.py` | **The Brain.** Manages `state.json` and `results.json`, ensuring data flows correctly between modules. |
| `utils.py` | Helper functions for system cleanup and miscellaneous tasks. |

---

## 🔬 How It Works

### 1. The Workflow Manager (`workflow.py`)
Instead of passing variables manually between notebook cells, the `WorkflowManager` saves the current state of the experiment to `state.json`.
* *Example:* The **Builder** saves the PDB path to the state. The **Simulator** reads that path from the state automatically.

### 2. Simulation Details
* **Force Field:** Generic OPLS (via OpenMM/RDKit generation).
* **Integrator:** Langevin Integrator (Simulating thermal bath collisions).
* **Conditions:** Vacuum simulation (No explicit solvent).

### 3. Visualization
The viewer (`visualiser.py`) uses **py3Dmol** but is enhanced with:
* **Custom Styling:** Standard CPK coloring with optimized sphere/stick ratios.
* **HUD:** A Heads-Up Display showing the current simulation step.
* **Frame Inspector:** A widget to scrub through individual frames for detailed inspection.

---

## 📊 Outputs

The simulation generates a folder named after your polymer (e.g., `PolyStyrene/`). Inside, you will find:

* `trajectory.dcd`: The raw motion data of the atoms.
* `polymer.pdb`: The starting structure.
* `results.json`: A summary of the calculated metrics (e.g., Average $R_g$).
* `workflow.log`: A text log of the operations performed.

---

## 🤝 Contributing

Contributions are welcome! Please fork the repository and submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
