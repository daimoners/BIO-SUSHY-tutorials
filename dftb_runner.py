import os
import subprocess
import shutil
import re
import sys
import time
import resource
from ase import io

def set_stack_limit():
    """Forces unlimited stack memory to prevent Fortran freezes."""
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    except ValueError:
        pass

def format_time(seconds):
    """Converts seconds to mm:ss format."""
    if seconds < 0: return "--:--"
    m, s = divmod(int(seconds), 60)
    return f"{m:02d}:{s:02d}"

def run_xtb_simulation(input_pdb, output_prefix, temp_k=300, steps=1000):
    print(f"--- ⚛️ Starting GFN2-xTB (Table Mode) ---")
    
    # 1. System Setup
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

    # 2. Prepare Input
    atoms = io.read(input_pdb)
    io.write("xtb_input.xyz", atoms)
    
    env = os.environ.copy()
    env["OMP_STACKSIZE"] = "4G"
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    
    # 3. Print Table Header
    print(f"   • Running {steps} steps...\n")
    
    # CSV-style Header matching your request
    # Step | Energy (eV) | Temp (K) | Speed (steps/s) | ETA
    header = f"{'Step':>8} | {'Energy (eV)':>14} | {'Temp (K)':>10} | {'Speed (stp/s)':>14} | {'ETA':>8}"
    print(header)
    print("-" * len(header))
    
    # 4. Launch MD
    md_cmd = [xtb_bin, "xtb_input.xyz", "--omd", "--gfn", "2", "--steps", str(steps), "--T", str(temp_k)]
    
    with open("xtb_md.log", "w") as log_f:
        process = subprocess.Popen(md_cmd, stdout=log_f, stderr=subprocess.STDOUT, env=env)

    # 5. Monitor Loop (Reading File)
    start_time = time.time()
    last_read_line = 0
    
    # Regex to find data lines: "   10   0.010   -45.23..."
    # We look for lines starting with an integer
    data_pattern = re.compile(r"^\s*(\d+)\s+([\d\.]+)\s+([\-\d\.]+)\s+[\d\.]+\s+([\d\.]+)")
    
    try:
        while process.poll() is None:
            time.sleep(1.5) # Update every 1.5 seconds
            
            if os.path.exists("xtb_md.log"):
                with open("xtb_md.log", "r") as f:
                    lines = f.readlines()
                    
                    # Process only new lines
                    if len(lines) > last_read_line:
                        new_lines = lines[last_read_line:]
                        last_read_line = len(lines)
                        
                        for line in new_lines:
                            # Try to parse data line
                            # xTB Format: Cycle | Time | Energy | Grad | Temp
                            parts = line.split()
                            if len(parts) >= 5 and parts[0].isdigit():
                                step = int(parts[0])
                                energy_hartree = float(parts[2])
                                temp = float(parts[4])
                                
                                # Calculations
                                elapsed = time.time() - start_time
                                speed = step / elapsed if elapsed > 0 else 0
                                remaining_steps = steps - step
                                eta_seconds = remaining_steps / speed if speed > 0 else 0
                                
                                # Convert Energy to eV (1 Eh = 27.2114 eV)
                                energy_ev = energy_hartree * 27.2114
                                
                                # PRINT ROW (Appended)
                                print(f"{step:>8} | {energy_ev:>14.4f} | {temp:>10.1f} | {speed:>14.2f} | {format_time(eta_seconds):>8}")
                                sys.stdout.flush()

    except KeyboardInterrupt:
        process.kill()
        print("\n🛑 Stopped by user.")
        return {}, None

    print("-" * len(header))
    print("   • MD Finished.")

    # 6. Check Outputs
    if not os.path.exists("xtbopt.xyz"):
        print("❌ Error: MD failed. Last 10 lines:")
        os.system("tail -n 10 xtb_md.log")
        return {}, None

    # 7. Final Properties
    print("   • Calculating final electronic properties...")
    with open("xtb_sp.log", "w") as f:
         subprocess.run([xtb_bin, "xtbopt.xyz", "--sp", "--gfn", "2"], stdout=f, env=env)
    
    # 8. Parse Results
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
