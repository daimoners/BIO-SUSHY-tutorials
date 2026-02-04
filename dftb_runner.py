import os
import shutil
import subprocess
import sys
import re
import numpy as np
from ase import io
from IPython.display import clear_output

def run_dftb_optimization(input_pdb, output_dir, solvent="water"):
    """
    Runs GFN2-xTB with IN-PLACE status updates.
    Parses 'CYCLE N' and '* total energy' from verbose logs.
    """
    # 1. Setup Directories
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    xtb_bin = shutil.which("xtb")
    if not xtb_bin:
        possible_loc = os.path.abspath("xtb-6.6.1/bin/xtb")
        if os.path.exists(possible_loc): xtb_bin = possible_loc
    if not xtb_bin: return {"status": "failed", "error": "xTB binary not found."}

    results = {"status": "failed"}
    original_dir = os.getcwd()

    try:
        os.chdir(output_dir)
        
        # 2. Clean Input
        try:
            atoms = io.read(input_pdb)
            io.write("input.xyz", atoms)
        except Exception as e:
            return {"status": "failed", "error": f"PDB->XYZ Conversion failed: {e}"}
        
        # 3. Run Optimization (Live Updates)
        cmd_opt = [xtb_bin, "input.xyz", "--opt", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none": cmd_opt.extend(["--alpb", solvent])
            
        print("⏳ Initializing xTB...")
        
        current_step = 0
        current_energy = "Calculating..."
        
        with open("dftb_opt.log", "w") as log_file:
            process = subprocess.Popen(cmd_opt, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            
            # Read line by line
            for line in process.stdout:
                log_file.write(line)
                
                # --- Parse Step Number ---
                # Pattern: ...... CYCLE   33 ......
                if "CYCLE" in line:
                    match = re.search(r"CYCLE\s+(\d+)", line)
                    if match:
                        current_step = match.group(1)
                        # Update Display immediately
                        clear_output(wait=True)
                        print(f"🚀 DFTB Optimization Running...\n   ► Step:   {current_step}\n   ► Energy: {current_energy}")

                # --- Parse Energy ---
                # Pattern: * total energy  :  -128.2708977 Eh
                if "* total energy" in line:
                    match = re.search(r"total energy\s+:\s+([\-\d\.]+)", line)
                    if match:
                        current_energy = f"{match.group(1)} Eh"
                        # Update Display with new energy
                        clear_output(wait=True)
                        print(f"🚀 DFTB Optimization Running...\n   ► Step:   {current_step}\n   ► Energy: {current_energy}")
            
            process.wait()

        # Clear the "Running" message so only results remain
        clear_output(wait=True)

        # Handle Output File
        if not os.path.exists("xtbopt.xyz"):
            if os.path.exists("xtbopt.log"): shutil.copy("xtbopt.log", "xtbopt.xyz")
            else:
                files = [f for f in os.listdir() if f.endswith('.xyz') and 'input' not in f]
                if files: shutil.copy(files[0], "xtbopt.xyz")

        if not os.path.exists("xtbopt.xyz"):
            raise RuntimeError("Optimization finished but 'xtbopt.xyz' not found.")

        # 4. Single Point (Silent)
        cmd_sp = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none": cmd_sp.extend(["--alpb", solvent])

        with open("properties.out", "w") as f:
            subprocess.run(cmd_sp, stdout=f, stderr=subprocess.STDOUT, text=True)

        # 5. Parse Data
        traj_energies = parse_optimization_log("dftb_opt.log")
        props = parse_properties_log("properties.out")
        
        final_energy = props["energy"]
        if final_energy == 0.0 and traj_energies:
            final_energy = traj_energies[-1]

        results = {
            "status": "success",
            "geometry_path": os.path.abspath("xtbopt.xyz"),
            "trajectory_energies": traj_energies,
            "total_energy": final_energy,
            "gap": props["gap"],
            "dipole": props["dipole"]
        }

    except Exception as e:
        results["error"] = str(e)
    finally:
        os.chdir(original_dir)
        
    return results

def parse_optimization_log(filename):
    """Parses energies from the verbose log format."""
    energies = []
    if not os.path.exists(filename): return energies
    
    with open(filename, 'r') as f:
        for line in f:
            # Pattern: * total energy  :  -128.2708977 Eh
            if "* total energy" in line:
                match = re.search(r"total energy\s+:\s+([\-\d\.]+)", line)
                if match:
                    try: energies.append(float(match.group(1)))
                    except: pass
    return energies

def parse_properties_log(filename):
    data = {"gap": 0.0, "dipole": 0.0, "energy": 0.0}
    if not os.path.exists(filename): return data
    with open(filename, 'r') as f:
        content = f.read()
    
    e_match = re.search(r"TOTAL ENERGY\s+([\-\d\.]+)\s+Eh", content)
    if e_match: data["energy"] = float(e_match.group(1))

    g_match = re.search(r"HOMO-LUMO GAP\s+([\d\.]+)\s+eV", content)
    if g_match: data["gap"] = float(g_match.group(1))

    lines = content.splitlines()
    for i, line in enumerate(lines):
        if "molecular dipole:" in line:
            for offset in range(1, 6):
                if i+offset < len(lines) and "full:" in lines[i+offset]:
                    try: data["dipole"] = float(lines[i+offset].split()[-1])
                    except: pass
                    break
    return data
