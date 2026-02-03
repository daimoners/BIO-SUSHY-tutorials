import os
import shutil
import subprocess
import numpy as np
from ase import io

def run_dftb_optimization(input_pdb, output_dir, solvent="water"):
    """
    Runs GFN2-xTB geometry optimization.
    - Robustly handles missing 'xtbopt.xyz' by checking for 'xtbopt.log'.
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
        
        # 4. Run Optimization
        cmd_opt = [xtb_bin, "input.xyz", "--opt", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none":
            cmd_opt.extend(["--alpb", solvent])
            
        # Capture output
        proc = subprocess.run(cmd_opt, capture_output=True, text=True)
        
        # Write log
        with open("dftb_opt.log", "w") as f:
            f.write(proc.stdout)
            f.write(proc.stderr)

        # --- FIX: Handle Missing .xyz File ---
        if not os.path.exists("xtbopt.xyz"):
            # Fallback 1: 'xtbopt.log' often contains the XYZ coordinates
            if os.path.exists("xtbopt.log"):
                shutil.copy("xtbopt.log", "xtbopt.xyz")
            # Fallback 2: Check for ANY .xyz file created recently
            else:
                files = [f for f in os.listdir() if f.endswith('.xyz') and 'input' not in f]
                if files:
                    shutil.copy(files[0], "xtbopt.xyz")

        # Double Check
        if not os.path.exists("xtbopt.xyz"):
            tail = "\n".join(proc.stdout.splitlines()[-20:])
            # List files to help debug
            file_list = "\n".join(os.listdir())
            raise RuntimeError(f"xTB finished but no geometry file found.\nFiles in folder:\n{file_list}\n\nxTB Log Tail:\n{tail}")

        # 5. Property Calculation
        cmd_sp = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2", "--chrg", "0"]
        if solvent and solvent != "none":
            cmd_sp.extend(["--alpb", solvent])

        with open("properties.out", "w") as f:
            subprocess.run(cmd_sp, stdout=f, stderr=subprocess.STDOUT, text=True)

        # 6. Parse Data
        gap, dipole = parse_properties("properties.out")
        charges = parse_charges("charges")
        polarity = np.std(charges) if charges else 0.0
        
        results = {
            "status": "success",
            "geometry_path": os.path.abspath("xtbopt.xyz"),
            "log_path": os.path.abspath("properties.out"),
            "band_gap": gap,
            "dipole": dipole,
            "polarity_index": polarity,
            "charges": charges
        }

    except Exception as e:
        results["error"] = str(e)
    finally:
        os.chdir(original_dir)
        
    return results

def parse_properties(filename):
    gap, dipole = 0.0, 0.0
    if not os.path.exists(filename): return gap, dipole
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    for i, line in enumerate(lines):
        if "HOMO-LUMO GAP" in line:
            try: gap = float(line.split()[-2])
            except: pass
        if "molecular dipole:" in line:
            for offset in range(1, 6):
                if i+offset < len(lines) and "full:" in lines[i+offset]:
                    try: dipole = float(lines[i+offset].split()[-1])
                    except: pass
                    break
    return gap, dipole

def parse_charges(filename):
    if os.path.exists(filename):
        with open(filename) as f:
            try: return [float(l.strip()) for l in f]
            except: return []
    return []
