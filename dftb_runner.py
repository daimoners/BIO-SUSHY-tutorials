import os
import shutil
import subprocess
import sys
import numpy as np
from ase import io

def run_dftb_optimization(input_pdb, output_dir, solvent="water"):
    """
    Runs GFN2-xTB geometry optimization with REAL-TIME output streaming.
    """
    # 1. Setup Directories
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    # 2. Locate xTB
    xtb_bin = shutil.which("xtb")
    if not xtb_bin:
        possible_loc = os.path.abspath("xtb-6.6.1/bin/xtb")
        if os.path.exists(possible_loc):
            xtb_bin = possible_loc
    
    if not xtb_bin:
        return {"status": "failed", "error": "xTB binary not found."}

    results = {"status": "failed"}
    original_dir = os.getcwd()

    try:
        os.chdir(output_dir)
        
        # 3. Clean Input: PDB -> XYZ
        try:
            atoms = io.read(input_pdb)
            io.write("input.xyz", atoms)
        except Exception as e:
            return {"status": "failed", "error": f"Failed to convert PDB to XYZ: {e}"}
        
        # 4. Run Optimization (STREAMED)
        print(f"   ► Running xTB Optimization in: {output_dir}")
        cmd_opt = [xtb_bin, "input.xyz", "--opt", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none":
            cmd_opt.extend(["--alpb", solvent])
            
        # Use Popen to stream output line by line
        with open("dftb_opt.log", "w") as log_file:
            process = subprocess.Popen(cmd_opt, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            
            # Stream lines
            for line in process.stdout:
                sys.stdout.write(line) # Print to screen
                log_file.write(line)   # Save to file
            
            process.wait() # Wait for finish

        # Check for success
        if not os.path.exists("xtbopt.xyz"):
            if os.path.exists("xtbopt.log"):
                shutil.copy("xtbopt.log", "xtbopt.xyz")
            else:
                files = [f for f in os.listdir() if f.endswith('.xyz') and 'input' not in f]
                if files: shutil.copy(files[0], "xtbopt.xyz")

        if not os.path.exists("xtbopt.xyz"):
            raise RuntimeError("Optimization finished but no geometry file found.")

        # 5. Single Point (Silent, it's fast)
        cmd_sp = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none":
            cmd_sp.extend(["--alpb", solvent])

        with open("properties.out", "w") as f:
            subprocess.run(cmd_sp, stdout=f, stderr=subprocess.STDOUT, text=True)

        # 6. Parse Data
        metrics = parse_xtb_output("properties.out")
        
        results = {
            "status": "success",
            "geometry_path": os.path.abspath("xtbopt.xyz"),
            "total_energy": metrics.get("energy", 0.0),
            "gap": metrics.get("gap", 0.0),
            "dipole": metrics.get("dipole", 0.0)
        }

    except Exception as e:
        results["error"] = str(e)
    finally:
        os.chdir(original_dir)
        
    return results

def parse_xtb_output(filename):
    data = {"gap": 0.0, "dipole": 0.0, "energy": 0.0}
    if not os.path.exists(filename): return data
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        if "TOTAL ENERGY" in line and "Eh" in line:
            try: data["energy"] = float(line.split()[-2])
            except: pass
        if "HOMO-LUMO GAP" in line:
            try: data["gap"] = float(line.split()[-2])
            except: pass
        if "molecular dipole:" in line:
            for offset in range(1, 6):
                if i+offset < len(lines) and "full:" in lines[i+offset]:
                    try: data["dipole"] = float(lines[i+offset].split()[-1])
                    except: pass
                    break
    return data
