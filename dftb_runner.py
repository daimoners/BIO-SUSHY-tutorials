import os
import subprocess
import shutil
import re
from ase import io

def find_xtb_executable():
    """Locates the xtb binary in the system."""
    # Common locations in Conda/Colab
    candidates = [
        shutil.which("xtb"),                  # Standard PATH
        "/usr/local/bin/xtb",                 # CondaColab default
        "/opt/conda/bin/xtb",                 # Standard Conda
        "/usr/bin/xtb"
    ]
    for cmd in candidates:
        if cmd and os.path.exists(cmd):
            return cmd
    return None

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    """
    Runs Quantum MD using the xTB binary directly.
    """
    print(f"--- ⚛️ Starting GFN2-xTB (Robust Mode) ---")
    
    # 1. Locate Binary
    xtb_bin = find_xtb_executable()
    if not xtb_bin:
        print("❌ Error: Could not find 'xtb' executable.")
        print("   -> Did you run the 'Fix Environment' cell?")
        return {}, None
    print(f"   • Found xTB binary at: {xtb_bin}")

    # 2. Prepare Input
    if not os.path.exists(input_pdb):
        print(f"❌ Error: Input file not found: {input_pdb}")
        return {}, None

    atoms = io.read(input_pdb)
    input_xyz = "xtb_input.xyz"
    io.write(input_xyz, atoms)
    print(f"   • Converted PDB to XYZ: {input_xyz}")

    # 3. Run MD (Capture output to see errors)
    print(f"   • Running {steps} steps of MD...")
    
    # We set OMP_NUM_THREADS to avoid Colab crashes
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"

    md_cmd = [xtb_bin, input_xyz, "--omd", "--gfn", "2", "--steps", str(steps), "--T", str(temp_k)]
    
    with open("xtb_md.log", "w") as log_f:
        try:
            # We run without 'shell=True' for better safety with the explicit binary path
            process = subprocess.run(md_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env, check=True)
        except subprocess.CalledProcessError:
            print("❌ xTB MD crashed! Last 10 lines of log:")
            print("-" * 40)
            os.system("tail -n 10 xtb_md.log")
            print("-" * 40)
            return {}, None

    # 4. Check Outputs
    if not os.path.exists("xtbopt.xyz"):
        print("❌ Error: MD finished but 'xtbopt.xyz' was not created.")
        return {}, None

    # 5. Final Properties Calculation
    print(f"   • Calculating final electronic properties...")
    sp_cmd = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"]
    
    with open("xtb_sp.log", "w") as log_f:
        subprocess.run(sp_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env, check=True)

    # 6. Parse Results
    results = parse_xtb_log("xtb_sp.log")
    
    # Save PDB
    final_atoms = io.read("xtbopt.xyz")
    final_pdb = f"{output_prefix}_final.pdb"
    io.write(final_pdb, final_atoms)
    
    # Cleanup Trajectory
    if os.path.exists("xtb.trj"):
        shutil.move("xtb.trj", f"{output_prefix}_traj.xyz")

    print(f"   ✅ Success. Gap: {results.get('homo_lumo_gap_eV', 0):.3f} eV")
    return results, final_pdb

def parse_xtb_log(logfile):
    """Regex parser for properties."""
    results = {}
    try:
        with open(logfile, 'r') as f:
            content = f.read()
            
        gap_match = re.search(r"HOMO-LUMO gap\s+:\s+([\d\.]+)\s+eV", content)
        dip_match = re.search(r"molecular dipole:\s+([\d\.]+)\s+Debye", content)
        
        if gap_match: results["homo_lumo_gap_eV"] = float(gap_match.group(1))
        if dip_match: results["dipole_moment_debye"] = float(dip_match.group(1))
    except Exception as e:
        print(f"⚠️ Error parsing logs: {e}")
        
    return results
