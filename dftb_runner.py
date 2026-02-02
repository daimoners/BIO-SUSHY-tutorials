import os
import subprocess
import shutil
import re
import sys
import time
from ase import io

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    """
    Runs xTB with a robust 'File Watcher' to ensure output is visible.
    """
    print(f"--- ⚛️ Starting GFN2-xTB (Watcher Mode) ---")
    
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

    # 2. Prepare Input & Check Size
    atoms = io.read(input_pdb)
    num_atoms = len(atoms)
    print(f"   • System Size: {num_atoms} atoms")
    
    if num_atoms > 300:
        print("   ⚠️ WARNING: System is large for GFN2-xTB on CPU.")
        print("   ⚠️ This might take 1-2 seconds PER STEP (Total: ~15 mins).")
        
    input_xyz = "xtb_input.xyz"
    io.write(input_xyz, atoms)
    
    # 3. Launch MD in Background
    print(f"   • Launching MD ({steps} steps)...")
    
    md_cmd = [xtb_bin, input_xyz, "--omd", "--gfn", "2", "--steps", str(steps), "--T", str(temp_k)]
    
    # Force single thread for stability if it was hanging
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1" 
    env["MKL_NUM_THREADS"] = "1"

    # Open the log file for writing
    with open("xtb_md.log", "w") as log_f:
        # Start the process independent of Python's buffer
        process = subprocess.Popen(md_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env)

    # 4. Watch the Log File
    start_time = time.time()
    step_pattern = re.compile(r"cycle\s+(\d+)")
    last_step = 0
    
    try:
        while process.poll() is None:
            # Sleep briefly
            time.sleep(2.0)
            
            # Read the log file to see progress
            if os.path.exists("xtb_md.log"):
                with open("xtb_md.log", "r") as f:
                    lines = f.readlines()
                    
                    # If log is empty, it might be initializing
                    if not lines:
                        continue

                    # Parse the last few lines for "cycle"
                    found_step = False
                    for line in reversed(lines[-10:]): # Check last 10 lines
                        match = step_pattern.search(line)
                        if match:
                            current_step = int(match.group(1))
                            last_step = current_step
                            
                            # Calc time
                            elapsed = time.time() - start_time
                            speed = elapsed / current_step if current_step > 0 else 0
                            eta = (steps - current_step) * speed
                            
                            sys.stdout.write(f"\r     ⏳ Step {current_step}/{steps} | ETA: {int(eta)}s | Last Log: {line.strip()[:40]}...")
                            sys.stdout.flush()
                            found_step = True
                            break
                    
                    # If no step found yet, print initialization info
                    if not found_step and len(lines) > 0:
                        sys.stdout.write(f"\r     ⚙️ Initializing... (xTB is running setup)")
                        sys.stdout.flush()
                        
    except KeyboardInterrupt:
        print("\n🛑 User stopped simulation.")
        process.kill()
        return {}, None

    print("\n   • MD Finished.")

    # 5. Check if it actually worked
    if not os.path.exists("xtbopt.xyz"):
        print("❌ Error: MD finished but output is missing.")
        print("   -> Printing last 20 lines of log for debugging:")
        print("-" * 50)
        os.system("tail -n 20 xtb_md.log")
        print("-" * 50)
        return {}, None

    # 6. Final Properties (Single Point)
    print(f"   • Calculating properties...")
    sp_cmd = [xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"]
    with open("xtb_sp.log", "w") as log_f:
        subprocess.run(sp_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env)

    results = parse_xtb_log("xtb_sp.log")
    
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
