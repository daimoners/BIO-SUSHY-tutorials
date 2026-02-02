import os
import subprocess
import shutil
import re
import sys
import time
from ase import io

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    """
    Runs Quantum MD using xTB.
    Features: ROBUST progress parsing (fixes 'ANC' error) + Speedometer.
    """
    print(f"--- ⚛️ Starting GFN2-xTB (Manual Mode) ---")
    
    # 1. Setup Binary
    local_xtb = os.path.abspath("xtb-6.6.1/bin/xtb")
    if os.path.exists(local_xtb):
        xtb_bin = local_xtb
        os.environ["XTBPATH"] = os.path.abspath("xtb-6.6.1/share/xtb")
    else:
        xtb_bin = shutil.which("xtb")
        
    if not xtb_bin:
        print("❌ Error: Could not find 'xtb' binary.")
        return {}, None

    # 2. Prepare Input
    atoms = io.read(input_pdb)
    input_xyz = "xtb_input.xyz"
    io.write(input_xyz, atoms)
    
    # 3. Run MD with LIVE Regex Parsing
    print(f"   • Running {steps} steps of MD...")
    
    md_cmd = [xtb_bin, input_xyz, "--omd", "--gfn", "2", "--steps", str(steps), "--T", str(temp_k)]
    env = os.environ.copy() 
    
    start_time = time.time()

    with open("xtb_md.log", "w") as log_f:
        process = subprocess.Popen(
            md_cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            env=env,
            text=True,
            bufsize=1
        )
        
        # Regex to find "cycle" followed by a number
        # Example xTB line: "     cycle      13 time ..."
        step_pattern = re.compile(r"cycle\s+(\d+)")
        
        for line in process.stdout:
            log_f.write(line)
            
            # Check for step number safely
            match = step_pattern.search(line)
            if match:
                current_step = int(match.group(1))
                
                # Calculate speed
                elapsed = time.time() - start_time
                if current_step > 0:
                    sec_per_step = elapsed / current_step
                    remaining = (steps - current_step) * sec_per_step
                    
                    # Print progress bar
                    sys.stdout.write(f"\r     ⏳ Step {current_step} / {steps} | Est. Remaining: {int(remaining)}s ")
                    sys.stdout.flush()
        
        process.wait()
        print("\n   • MD Finished.")

    # 4. Final Properties (Single Point)
    if not os.path.exists("xtbopt.xyz"):
        print("❌ Error: MD crashed. Check 'xtb_md.log'.")
        return {}, None

    print(f"   • Calculating electronic properties (Final SP)...")
    sp_cmd = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"]
    
    with open("xtb_sp.log", "w") as log_f:
        subprocess.run(sp_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env, check=True)

    # 5. Parse Results
    results = parse_xtb_log("xtb_sp.log")
    
    # Save files
    final_atoms = io.read("xtbopt.xyz")
    final_pdb = f"{output_prefix}_final.pdb"
    io.write(final_pdb, final_atoms)
    
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
