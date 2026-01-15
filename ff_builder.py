import os
import shutil
import warnings
import parmed as pmd
import MDAnalysis as mda
from foyer import Forcefield

# ==========================================
#        HELPER FUNCTIONS
# ==========================================

def sanitize_xml(xml_path):
    # Fixes XML attributes that crash the new OpenMM parser
    with open(xml_path, 'r') as f: lines = f.readlines()
    new_path = xml_path.replace(".xml", "_fixed.xml")
    with open(new_path, 'w') as f:
        for line in lines:
            if "combining_rule" in line:
                line = line.replace('combining_rule="geometric"', '').replace("combining_rule='geometric'", '')
            f.write(line)
    return new_path

def read_top_file(top_file):
    with open(top_file, "r") as top:
        return top.readlines()

def assign_foyer_charges(topology, num_atoms):
    """ Extract charges from the generated topology """
    line = 0
    while line < len(topology) and "[ atoms ]" not in topology[line]:
        line += 1
    while line < len(topology) and "opls_" not in topology[line]:
        line += 1
    charges = []
    for _ in range(num_atoms):
        if line >= len(topology): break
        parts = topology[line].split()
        if len(parts) > 6:
            charges.append(parts[6])
        line += 1
    return charges

def neutralize_charges(charges):
    """ Smear residual charge to ensure net zero """
    tot_charge = sum([float(x) for x in charges])
    if abs(tot_charge) < 1e-5: return charges
    
    dq_c = tot_charge / len(charges)
    return [str(float(charge) - dq_c) for charge in charges]

def write_neutralized_topology(top_lines, out_path, charges, num_atoms):
    with open(out_path, 'w') as f:
        line = 0
        while line < len(top_lines) and "opls_" not in top_lines[line]:
            f.write(top_lines[line])
            line += 1
        
        atom_count = 0
        while atom_count < num_atoms and line < len(top_lines):
            parts = top_lines[line].split()
            if len(parts) > 6:
                parts[6] = f"{float(charges[atom_count]):.5f}"
                f.write("\t".join(parts) + "\n")
                atom_count += 1
            else:
                f.write(top_lines[line])
            line += 1
            
        while line < len(top_lines):
            f.write(top_lines[line])
            line += 1

def inspect_system(gro_file, top_file):
    """
    Prints a structured summary of the generated files 
    so the user can verify the parameterization.
    """
    print("\n" + "="*60)
    print(f"👀  INSPECTING GENERATED FORCE FIELD")
    print("="*60)

    # 1. Show Structure Head (GRO)
    print(f"\n📂 Structure File ({os.path.basename(gro_file)}):")
    with open(gro_file, 'r') as f:
        for i in range(5): 
            print(f"   {f.readline().strip()}")
    print("   ... (coordinates continue) ...")

    # 2. Show Topology Details (TOP)
    print(f"\n📜 Topology File ({os.path.basename(top_file)}):")
    
    with open(top_file, 'r') as f:
        lines = f.readlines()

    print_section = False
    printed_count = 0
    sections_to_show = ["[ atomtypes ]", "[ atoms ]", "[ bonds ]"]

    for line in lines:
        line = line.strip()
        
        # Detect Header
        if line.startswith("[") and line.endswith("]"):
            if any(sec in line for sec in sections_to_show):
                print(f"\n👉 {line}")
                print_section = True
                printed_count = 0
                continue
            else:
                print_section = False

        # Print Content
        if print_section and line and not line.startswith(";"):
            if printed_count < 6:
                print(f"   {line}")
                printed_count += 1
            elif printed_count == 6:
                print("   ... (truncated) ...")
                printed_count += 1
    
    print("-" * 60)
    print("✅ System ready for Simulation.")

# ==========================================
#        MAIN BUILDER FUNCTION
# ==========================================

# FIX: Added output_name argument back here
def build_opls_system(pdb_file, output_name="system_opls"):
    print(f"⚙️ OPLS BUILDER: Parameterizing {pdb_file}...")

    # 1. Locate and Fix XML
    base_dir = os.path.dirname(os.path.abspath(__file__))
    xml_path = os.path.join(base_dir, "oplsaa.xml")
    if not os.path.exists(xml_path):
        xml_path = "oplsaa.xml" # Fallback

    clean_xml = sanitize_xml(xml_path)

    # 2. Load Molecule
    mol = pmd.load_file(pdb_file)

    # CRITICAL FIX: Add virtual simulation box (Foyer requires this)
    if mol.box is None:
        print("   -> Adding virtual simulation box (10 nm)...")
        mol.box = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]

    # 3. Apply Force Field
    print("   -> Atom-typing with Foyer...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ff = Forcefield(forcefield_files=clean_xml)
        typed_structure = ff.apply(mol)

    # 4. Save Output
    out_dir = f"{output_name}"
    os.makedirs(out_dir, exist_ok=True)

    gro_path = os.path.join(out_dir, "conf.gro")
    top_path = os.path.join(out_dir, "topol.top")

    typed_structure.save(gro_path, overwrite=True)
    typed_structure.save(top_path, overwrite=True)

    print(f"✅ Success! OPLS System saved to: {out_dir}")
    return gro_path, top_path
