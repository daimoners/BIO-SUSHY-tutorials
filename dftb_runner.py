import os
import subprocess
import shutil
import re
import sys
import time
import resource # <--- NEW: To control system limits
from ase import io

def set_stack_limit():
    """Forces the system to allow unlimited stack memory (fixes Fortran hangs)."""
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except ValueError:
        pass # Sometimes Colab forbids this, but we try anyway

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    print(f"--- ⚛️ Starting GFN2-xTB (Debug Mode) ---")
    
    # 1. System Setup
    set_stack_limit() # Apply memory fix
    
    local_xtb = os.path.abspath("xtb-6.6.1/bin/xtb")
    if os.path.exists(local_xtb):
        xtb_bin = local_xtb
        os.environ["XTBPATH"] = os.path.abspath("xtb-6.6.1/share/xtb")
    else:
        print("❌ Error: xTB binary not found.")
        return {}, None

    # 2. Prepare Input
    atoms = io.read(input_pdb)
    io.write("xtb_input.xyz", atoms)
    
    # 3. Environment Variables (CRITICAL FOR STABILITY)
    env = os.environ.copy()
    env["OMP_STACKSIZE"] = "4G"  # Give it 4GB of stack memory
    env["OMP_NUM_THREADS"] = "1" # Force serial execution (stops threading hangs)
    env["MKL_NUM_THREADS"] = "1"
    
    print(f"   • Running {steps} steps (Safe Mode)...")
    
    # 4. Launch with Live Output
    md_cmd = [xtb_bin, "xtb_input.xyz", "--omd", "--gfn", "2", "--steps", str(steps), "--T", str(temp_k)]
    
    with open("xtb_md.log", "w") as log_f:
        process = subprocess.Popen(md_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env)

    # 5. Monitor Loop
    start_time = time.time()
    step_pattern = re.compile(r"cycle\s+(\d+)")
    
    try:
        while process.poll() is None:
            time.sleep(2.0)
            
            if os.path.exists("xtb_md.log"):
                with open("xtb_md.log", "r") as f:
                    lines = f.readlines()
                    if not lines: continue

                    # Check for progress
                    found_step = False
                    for line in reversed(lines[-20:]):
                        match = step_pattern.search(line)
                        if match:
                            step = int(match.group(1))
                            sys.stdout.write(f"\r     ⏳ Progress: Step {step}/{steps} ")
                            sys.stdout.flush()
                            found_step = True
                            break
                    
                    # IF NO PROGRESS: Print the raw log so we see why it's stuck
                    if not found_step:
                        last_line = lines[-1].strip()
                        if len(last_line) > 50: last_line = last_line[:50] + "..."
                        sys.stdout.write(f"\r     ⚙️ Status: {last_line}")
                        sys.stdout.flush()

    except KeyboardInterrupt:
        process.kill()
        print("\n🛑 Stopped by user.")
        return {}, None

    print("\n   • MD Finished.")
    
    # 6. Check for output
    if not os.path.exists("xtbopt.xyz"):
        print("❌ Error: MD failed. Last 10 lines:")
        os.system("tail -n 10 xtb_md.log")
        return {}, None

    # 7. Final Calc
    print("   • Calculating properties...")
    subprocess.run([xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    
    # 8. Parse
    results = {}
    if os.path.exists("xtbopt.log"): # Default name for SP calc might vary, checking main log
        # Usually SP writes to stdout. Let's run explicitly to a file.
        with open("xtb_sp.log", "w") as f:
             subprocess.run([xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"], stdout=f, env=env)
        
        with open("xtb_sp.log", "r") as f: content = f.read()
        gap = re.search(r"HOMO-LUMO gap\s+:\s+([\d\.]+)\s+eV", content)
        dip = re.search(r"molecular dipole:\s+([\d\.]+)\s+Debye", content)
        if gap: results["homo_lumo_gap_eV"] = float(gap.group(1))
        if dip: results["dipole_moment_debye"] = float(dip.group(1))

    final_pdb = f"{output_prefix}_final.pdb"
    io.write(final_pdb, io.read("xtbopt.xyz"))
    
    print(f"   ✅ Done. Gap: {results.get('homo_lumo_gap_eV', 0):.3f} eV")
    return results, final_pdb
