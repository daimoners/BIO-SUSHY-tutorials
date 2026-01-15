import os
import shutil
import parmed as pmd
import MDAnalysis as mda
from foyer import Forcefield

# ==========================================
#        HELPER FUNCTIONS
# ==========================================

def sanitize_xml(xml_path):
    """
    Fixes compatibility issues between older XMLs and newer Foyer/OpenMM versions.
    Specifically removes the 'combining_rule' attribute from the <ForceField> tag.
    """
    temp_path = xml_path.replace(".xml", "_fixed.xml")
    
    with open(xml_path, 'r') as f:
        lines = f.readlines()
    
    with open(temp_path, 'w') as f:
        for line in lines:
            if "<ForceField" in line and "combining_rule" in line:
                # Remove the offending attribute
                line = line.replace('combining_rule="geometric"', '')
                line = line.replace("combining_rule='geometric'", '')
            f.write(line)
            
    return temp_path

def read_top_file(top_file):
    with open(top_file, "r") as top:
        return top.readlines()

def assign_foyer_charges(topology, num_atoms):
    """ Extract charges from the [ atoms ] section of a Foyer-generated topology """
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
    
    print(f"      ! Neutralizing small net charge: {tot_charge:.4f} e")
    dq_c = tot_charge / len(charges)
    return [str(float(charge) - dq_c) for charge in charges]

def write_neutralized_topology(top_lines, out_path, charges, num_atoms):
    with open(out_path, 'w') as f:
        line = 0
        # Write headers
        while line < len(top_lines) and "opls_" not in top_lines[line]:
            f.write(top_lines[line])
            line += 1
        
        # Write Atoms with new charges
        atom_count = 0
        while atom_count < num_atoms and line < len(top_lines):
            parts = top_lines[line].split()
            if len(parts) > 6:
                # Format: nr type resnr residue atom cgnr charge mass
                # We replace the charge (index 6)
                parts[6] = f"{float(charges[atom_count]):.5f}"
                f.write("\t".join(parts) + "\n")
                atom_count += 1
            else:
                f.write(top_lines[line])
            line += 1
            
        # Write the rest
        while line < len(top_lines):
            f.write(top_lines[line])
            line += 1

# ==========================================
#        MAIN BUILDER FUNCTION
# ==========================================

def build_opls_system(pdb_file, resname="POL"):
    """
    Applies OPLS-AA Force Field using Foyer and writes GROMACS files.
    """
    # 1. Setup Paths
    base_dir = os.path.dirname(os.path.abspath(__file__))
    original_xml = os.path.join(base_dir, "oplsaa.xml")
    
    if not os.path.exists(original_xml):
        raise FileNotFoundError(f"Could not find oplsaa.xml at: {original_xml}")

    # 2. Sanitize XML (Fix for Colab/Newer Foyer)
    ff_xml = sanitize_xml(original_xml)

    # Create output directory
    job_name = os.path.basename(pdb_file).replace(".pdb", "")
    output_dir = os.path.abspath(f"{job_name}_OPLS")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⚙️ OPLS BUILDER: Processing {pdb_file}...")
    print(f"   -> Output Directory: {output_dir}")

    # 3. Load Structure into ParmEd
    try:
        mol = pmd.load_file(pdb_file)
    except Exception as e:
        raise ValueError(f"Could not load PDB: {e}")

    # 4. Apply Force Field (Foyer)
    print("   -> Atom-typing with Foyer (this may take a moment)...")
    try:
        # Load the sanitized XML
        ff = Forcefield(forcefield_files=ff_xml)
        typed_structure = ff.apply(mol)
    except Exception as e:
        print(f"\n❌ Foyer failed to atom-type the molecule.")
        raise e

    # 5. Save Intermediate GROMACS files
    top_file = os.path.join(output_dir, "foyer_out.top")
    gro_file = os.path.join(output_dir, "conf.gro")
    
    typed_structure.save(top_file, overwrite=True)
    typed_structure.save(gro_file, overwrite=True)
    
    # 6. Post-Processing (Neutralization & Cleanup)
    print("   -> Post-processing topology (Charge Neutralization)...")
    
    top_lines = read_top_file(top_file)
    u = mda.Universe(gro_file)
    n_atoms = len(u.atoms)
    
    raw_charges = assign_foyer_charges(top_lines, n_atoms)
    clean_charges = neutralize_charges(raw_charges)
    
    # Write Final ITP/TOP
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

    # Cleanup temp XML
    if os.path.exists(ff_xml): os.remove(ff_xml)

    print(f"✅ OPLS System Ready!")
    return gro_file, final_itp, final_top
