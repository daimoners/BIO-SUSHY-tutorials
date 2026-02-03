%%writefile dftb_runner.py
import os
import shutil
import subprocess
import numpy as np

def run_dftb_optimization(input_pdb, output_dir, solvent="water"):
    """
    Runs GFN2-xTB geometry optimization and property calculation.
    Returns: dictionary of results and paths.
    """
    # 1. Setup Directories
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    
    # 2. Locate xTB
    xtb_bin = shutil.which("xtb")
    if not xtb_bin:
        # Fallback search
        possible_loc = os.path.abspath("xtb-6.6.1/bin/xtb")
        if os.path.exists(possible_loc):
            xtb_bin = possible_loc
    
    if not xtb_bin:
        raise FileNotFoundError("xTB binary not found. Please install xTB.")

    results = {"status": "failed"}
    original_dir = os.getcwd()

    try:
        os.chdir(output_dir)
        
        # 3. Geometry Optimization
        shutil.copy(input_pdb, "input.pdb")
        cmd_opt = [xtb_bin, "input.pdb", "--opt", "--gfn", "2"]
        if solvent and solvent != "none":
            cmd_opt.extend(["--alpb", solvent])
            
        subprocess.run(cmd_opt, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        if not os.path.exists("xtbopt.xyz"):
            raise RuntimeError("Geometry optimization failed.")

        # 4. Property Calculation (Single Point)
        cmd_sp = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"]
        if solvent and solvent != "none":
            cmd_sp.extend(["--alpb", solvent])

        with open("properties.out", "w") as f:
            subprocess.run(cmd_sp, stdout=f, stderr=subprocess.STDOUT, text=True)

        # 5. Parse Data
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
            return [float(l.strip()) for l in f]
    return []
