%%writefile dftb_runner.py
import os
import shutil
import subprocess
import sys
import re
import numpy as np
from ase import io

def run_dftb_optimization(input_pdb, output_dir, solvent="water"):
    """
    Runs GFN2-xTB geometry optimization.
    - Streams output.
    - Captures Energy Trajectory for plotting.
    - Fixes "0.00 eV" parsing bug.
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
        
        # 3. Run Optimization (Streaming)
        print(f"   ► Running xTB Optimization in: {output_dir}")
        cmd_opt = [xtb_bin, "input.xyz", "--opt", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none": cmd_opt.extend(["--alpb", solvent])
            
        with open("dftb_opt.log", "w") as log_file:
            process = subprocess.Popen(cmd_opt, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                sys.stdout.write(line)
                log_file.write(line)
            process.wait()

        # Handle Output File (Rename log if needed)
        if not os.path.exists("xtbopt.xyz"):
            if os.path.exists("xtbopt.log"): shutil.copy("xtbopt.log", "xtbopt.xyz")
            else:
                files = [f for f in os.listdir() if f.endswith('.xyz') and 'input' not in f]
                if files: shutil.copy(files[0], "xtbopt.xyz")

        if not os.path.exists("xtbopt.xyz"):
            raise RuntimeError("Optimization finished but 'xtbopt.xyz' not found.")

        # 4. Run Single Point (Reliable Properties)
        cmd_sp = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none": cmd_sp.extend(["--alpb", solvent])

        with open("properties.out", "w") as f:
            subprocess.run(cmd_sp, stdout=f, stderr=subprocess.STDOUT, text=True)

        # 5. Parse Data
        # A. Get Trajectory from Optimization Log
        traj_energies = parse_optimization_log("dftb_opt.log")
        
        # B. Get Final Properties from Single Point Log
        props = parse_properties_log("properties.out")
        
        # If final energy is missing in props, take the last trajectory point
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
    """Extracts energy vs cycle number from xTB optimization log."""
    energies = []
    if not os.path.exists(filename): return energies
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # xTB Geometry Optimization Table looks like:
    #  cycle      energy      ...
    #     1    -50.12345      ...
    start_reading = False
    for line in lines:
        if "cycle" in line and "energy" in line:
            start_reading = True
            continue
        
        if start_reading:
            parts = line.split()
            # Valid line: "   1   -50.123   ..." (Int then Float)
            if len(parts) >= 2 and parts[0].isdigit():
                try:
                    e = float(parts[1])
                    energies.append(e)
                except: pass
            elif "Average" in line or "Geometry" in line:
                # Table ended
                start_reading = False
                
    return energies

def parse_properties_log(filename):
    data = {"gap": 0.0, "dipole": 0.0, "energy": 0.0}
    if not os.path.exists(filename): return data
    
    with open(filename, 'r') as f:
        content = f.read()
        lines = content.splitlines()

    # 1. Energy: | TOTAL ENERGY  -83.935... Eh |
    # Regex is safer than split()[-2]
    e_match = re.search(r"TOTAL ENERGY\s+([\-\d\.]+)\s+Eh", content)
    if e_match: data["energy"] = float(e_match.group(1))

    # 2. Gap: | HOMO-LUMO GAP   5.489... eV |
    g_match = re.search(r"HOMO-LUMO GAP\s+([\d\.]+)\s+eV", content)
    if g_match: data["gap"] = float(g_match.group(1))

    # 3. Dipole (Table Search)
    for i, line in enumerate(lines):
        if "molecular dipole:" in line:
            for offset in range(1, 6):
                if i+offset < len(lines):
                    subline = lines[i+offset]
                    if "full:" in subline:
                        try: data["dipole"] = float(subline.split()[-1])
                        except: pass
                        break
    return data
