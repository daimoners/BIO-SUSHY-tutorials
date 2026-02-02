import os
import subprocess
import shutil
import re
import sys
import time
import resource
from ase import io

def set_stack_limit():
    """Forces unlimited stack memory."""
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except ValueError:
        pass

def format_time(seconds):
    """Converts seconds to mm:ss."""
    if seconds < 0: return "--:--"
    m, s = divmod(int(seconds), 60)
    return f"{m:02d}:{s:02d}"

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    print(f"--- ⚛️ Starting GFN2-xTB (Table Mode) ---")
    
    # 1. Setup
    set_stack_limit()
    local_xtb = os.path.abspath("xtb-6.6.1/bin/xtb")
    
    if os.path.exists(local_xtb):
        xtb_bin = local_xtb
        os.environ["XTBPATH"] = os.path.abspath("xtb-6.6.1/share/xtb")
    else:
        xtb_bin = shutil.which("xtb")
        if not xtb_bin:
            print("❌ Error: xTB binary not found.")
            return {}, None

    # 2. Input
    atoms = io.read(input_pdb)
    io.write("xtb_input.xyz", atoms)
    
    env = os.environ.copy()
    env["OMP_STACKSIZE"] = "4G"
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    
    # 3. Header
    print(f"   • Running {steps} steps...\n")
    header = f"{'Step':>8} | {'Energy (eV)':>14} | {'Temp (K)':>10} | {'Speed (stp/s)':>14} | {'ETA':>8}"
    print(header)
    print("-" * len(header))
    
    # 4. Launch Process
    md_cmd = [xtb_bin, "xtb_input.xyz", "--omd", "--gfn", "2", "--steps", str(steps), "--T", str(temp_k)]
    
    with open("xtb_md.log", "w") as log_f:
        process = subprocess.Popen(md_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env)

    # 5. Smart Watcher Loop
    start_time = time.time()
    last_read_line = 0
    
    try:
        while process.poll() is None:
            time.sleep(1.0)
            
            if os.path.exists("xtb_md.log"):
                with open("xtb_md.log", "r") as f:
                    lines = f.readlines()
                    
                    if len(lines) > last_read_line:
                        new_lines = lines[last_read_line:]
                        last_read_line = len(lines)
                        
                        for line in new_lines:
                            parts = line.split()
                            
                            # SAFETY CHECK: Only process lines that are purely numeric data
                            if len(parts) >= 5 and parts[0].isdigit():
                                try:
                                    # Try to convert columns to numbers. 
                                    # If any part is text (like 'C'), this fails and goes to 'except'
                                    step = int(parts[0])
                                    time_val = float(parts[1]) # Verify this is a number
                                    energy_hartree = float(parts[2])
                                    temp = float(parts[4])

                                    # Calc Stats
                                    elapsed = time.time() - start_time
                                    speed = step / elapsed if elapsed > 0 else 0
                                    eta = (steps - step) / speed if speed > 0 else 0
                                    energy_ev = energy_hartree * 27.2114
                                    
                                    # Only print every 10th step (or slow steps) to avoid spamming
                                    # But always print the first few
                                    if step % 10 == 0 or step < 5:
                                        print(f"{step:>8} | {energy_ev:>14.4f} | {temp:>10.1f} | {speed:>14.2f} | {format_time(eta):>8}")
                                        sys.stdout.flush()
                                        
                                except ValueError:
                                    # This was an atom line (e.g. "1 C ..."), ignore it.
                                    continue

    except KeyboardInterrupt:
        process.kill()
        print("\n🛑 Stopped by user.")
        return {}, None

    print("-" * len(header))
    print("   • MD Finished.")

    # 6. Finalize
    if not os.path.exists("xtbopt.xyz"):
        print("❌ Error: MD output missing.")
        return {}, None

    print("   • Calculating final properties...")
    with open("xtb_sp.log", "w") as f:
         subprocess.run([xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"], stdout=f, env=env)
    
    results = {}
    with open("xtb_sp.log", "r") as f: content = f.read()
    gap = re.search(r"HOMO-LUMO gap\s+:\s+([\d\.]+)\s+eV", content)
    dip = re.search(r"molecular dipole:\s+([\d\.]+)\s+Debye", content)
    
    if gap: results["homo_lumo_gap_eV"] = float(gap.group(1))
    if dip: results["dipole_moment_debye"] = float(dip.group(1))

    final_pdb = f"{output_prefix}_final.pdb"
    io.write(final_pdb, io.read("xtbopt.xyz"))
    
    print(f"   ✅ Done. Gap: {results.get('homo_lumo_gap_eV', 0):.3f} eV")
    return results, final_pdb
