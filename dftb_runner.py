import os
import subprocess
import shutil
import re
from ase import io

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    """
    Runs Quantum MD using the LOCALLY DOWNLOADED xTB binary.
    """
    print(f"--- ⚛️ Starting GFN2-xTB (Manual Mode) ---")
    
    # 1. FIND THE BINARY (The one we just downloaded)
    # We look in the current directory first
    local_xtb = os.path.abspath("xtb-6.6.1/bin/xtb")
    
    if os.path.exists(local_xtb):
        xtb_bin = local_xtb
        # We must set the XTBPATH so it finds its parameter files
        os.environ["XTBPATH"] = os.path.abspath("xtb-6.6.1/share/xtb")
    else:
        # Fallback to system path if local download is missing
        xtb_bin = shutil.which("xtb")
        
    if not xtb_bin:
        print("❌ Error: Could not find 'xtb' binary.")
        print("   -> Please run the 'Direct Install' cell above.")
        return {}, None
        
    print(f"   • Using Engine: {xtb_bin}")

    # 2. Prepare Input
    if not os.path.exists(input_pdb):
        print(f"❌ Error: Input file not found: {input_pdb}")
        return {}, None

    atoms = io.read(input_pdb)
    input_xyz = "xtb_input.xyz"
    io.write(input_xyz, atoms)
    
    # 3. Run MD
    print(f"   • Running {steps} steps of MD...")
    
    # Command: xtb input.xyz --omd --gfn 2 ...
    md_cmd = [xtb_bin, input_xyz, "--omd", "--gfn", "2", "--steps", str(steps), "--T", str(temp_k)]
    
    env = os.environ.copy() # Pass our modified env vars
    
    with open("xtb_md.log", "w") as log_f:
        try:
            subprocess.run(md_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env, check=True)
        except subprocess.CalledProcessError:
            print("❌ xTB MD crashed. Log tail:")
            os.system("tail -n 10 xtb_md.log")
            return {}, None

    # 4. Final Optimization (Single Point)
    if not os.path.exists("xtbopt.xyz"):
        print("❌ Error: No output structure (xtbopt.xyz).")
        return {}, None

    print(f"   • Calculating electronic properties...")
    sp_cmd = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"]
    
    with open("xtb_sp.log", "w") as log_f:
        subprocess.run(sp_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env, check=True)

    # 5. Parse & Save
    results = parse_xtb_log("xtb_sp.log")
    
    final_atoms = io.read("xtbopt.xyz")
    final_pdb = f"{output_prefix}_final.pdb"
    io.write(final_pdb, final_atoms)
    
    # Cleanup
    if os.path.exists("xtb.trj"):
        shutil.move("xtb.trj", f"{output_prefix}_traj.xyz")

    print(f"   ✅ Done. Gap: {results.get('homo_lumo_gap_eV', 0):.3f} eV")
    return results, final_pdb

def parse_xtb_log(logfile):
    results = {}
    try:
        with open(logfile, 'r') as f: content = f.read()
        gap = re.search(r"HOMO-LUMO gap\s+:\s+([\d\.]+)\s+eV", content)
        dip = re.search(r"molecular dipole:\s+([\d\.]+)\s+Debye", content)
        if gap: results["homo_lumo_gap_eV"] = float(gap.group(1))
        if dip: results["dipole_moment_debye"] = float(dip.group(1))
    except: pass
    return results
