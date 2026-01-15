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
    """
    Creates a temporary XML compatible with Foyer/OpenMM in Colab.
    Removes attributes that newer parsers reject.
    """
    temp_path = xml_path.replace(".xml", "_fixed.xml")
    
    with open(xml_path, 'r') as f:
        lines = f.readlines()
    
    with open(temp_path, 'w') as f:
        for line in lines:
            # Fix: Remove 'combining_rule' from <ForceField>
            if "<ForceField" in line and "combining_rule" in line:
                line = line.replace('combining_rule="geometric"', '')
                line = line.replace("combining_rule='geometric'", '')
            f.write(line)
            
    return temp_path

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

# ==========================================
#        MAIN BUILDER FUNCTION
# ==========================================

def build_opls_system(pdb_file, output_name="system_opls", resname="POL"):
    # 1. Setup Paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    original_xml = os.path.join(base_dir, "oplsaa.xml")
    
    if not os.path.exists(original_xml):
        # Fallback if XML is in current dir
        if os.path.exists("oplsaa.xml"):
             original_xml = "oplsaa.xml"
        else:
             raise FileNotFoundError(f"Missing oplsaa.xml at: {original_xml}")

    # 2. Fix XML Issues
    ff_xml = sanitize_xml(original_xml)

    # 3. Prepare Output
    job_name = os.path.basename(pdb_file).replace(".pdb", "")
    # Use the output_name argument if provided, or default to folder naming
    output_dir = os.path.abspath(output_name)
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⚙️ OPLS BUILDER: Processing {pdb_file}...")

    # 4. Load & Fix Structure
    try:
        mol = pmd.load_file(pdb_file)
        
        # CRITICAL FIX: Add Simulation Box
        # Foyer requires a periodic box, even for vacuum. 
        # We add a 10nm (100A) cubic box.
        if mol.box is None:
            print("   -> Adding virtual simulation box (10nm)...")
            mol.box = [100.0, 100.0, 100.0, 90.0, 90.0, 90.0]
            
    except Exception as e:
        raise ValueError(f"Could not load PDB: {e}")

    # 5. Apply Force Field (Foyer)
    print("   -> Atom-typing with Foyer...")
    try:
        # We capture warnings to keep output clean
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ff = Forcefield(forcefield_files=ff_xml)
            typed_structure = ff.apply(mol)
            
    except Exception as e:
        print(f"\n❌ Foyer Error: {e}")
        print("   (This usually means OpenMM version is > 8.0.0 or Box is missing)")
        raise e

    # 6. Save Outputs
    top_file = os.path.join(output_dir, "foyer_out.top")
    gro_file = os.path.join(output_dir, "conf.gro")
    
    typed_structure.save(top_file, overwrite=True)
    typed_structure.save(gro_file, overwrite=True)
    
    # 7. Post-Process (Neutralize)
    print("   -> Post-processing (Neutralization)...")
    top_lines = read_top_file(top_file)
    u = mda.Universe(gro_file)
    n_atoms = len(u.atoms)
    
    clean_charges = neutralize_charges(assign_foyer_charges(top_lines, n_atoms))
    
    final_itp = os.path.join(output_dir, f"{resname}.itp")
    write_neutralized_topology(top_lines, final_itp, clean_charges, n_atoms)
    
    final_top = os.path.join(output_dir, "topol.top")
    with open(final_top, 'w') as f:
        f.write(f"; OPLS-AA Topology for {job_name}\n")
        f.write(f'#include "{os.path.basename(original_xml)}"\n') 
        f.write(f'#include "{os.path.basename(final_itp)}"\n\n')
        f.write("[ system ]\n")
        f.write(f"{job_name}\n\n")
        f.write("[ molecules ]\n")
        f.write(f"Other 1\n")

    # Cleanup
    if os.path.exists(ff_xml): os.remove(ff_xml)

    print(f"✅ OPLS System Ready: {output_dir}")
    return gro_file, final_top
