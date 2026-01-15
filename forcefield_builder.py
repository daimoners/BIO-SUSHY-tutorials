import os
import shutil
import parmed as pmd
import MDAnalysis as mda
from foyer import Forcefield

# ==========================================
#        HELPER FUNCTIONS
# ==========================================

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
    """ Smear residual charge (common in capped polymers) to ensure net zero """
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
    ff_xml = os.path.join(base_dir, "oplsaa.xml")
    
    if not os.path.exists(ff_xml):
        raise FileNotFoundError(f"Could not find oplsaa.xml at: {ff_xml}")

    # Create output directory based on input name
    job_name = os.path.basename(pdb_file).replace(".pdb", "")
    output_dir = os.path.abspath(f"{job_name}_OPLS")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"⚙️ OPLS BUILDER: Processing {pdb_file}...")
    print(f"   -> Using Force Field: {os.path.basename(ff_xml)}")
    print(f"   -> Output Directory: {output_dir}")

    # 2. Load Structure into ParmEd
    try:
        mol = pmd.load_file(pdb_file)
        # Foyer requires a PDB with valid element info. 
        # Sometimes RDKit PDBs are messy, reloading via MDA or PMD helps.
    except Exception as e:
        raise ValueError(f"Could not load PDB: {e}")

    # 3. Apply Force Field (Foyer)
    print("   -> Atom-typing with Foyer (this may take a moment)...")
    try:
        ff = Forcefield(forcefield_files=ff_xml)
        typed_structure = ff.apply(mol)
    except Exception as e:
        print(f"\n❌ Foyer failed to atom-type the molecule.")
        print("   Common causes: Missing atom types in oplsaa.xml or bad geometry.")
        raise e

    # 4. Save Intermediate GROMACS files
    top_file = os.path.join(output_dir, "foyer_out.top")
    gro_file = os.path.join(output_dir, "conf.gro")
    
    typed_structure.save(top_file, overwrite=True)
    typed_structure.save(gro_file, overwrite=True)
    
    # 5. Post-Processing (Neutralization & Cleanup)
    print("   -> Post-processing topology (Charge Neutralization)...")
    
    # Read generated topology
    top_lines = read_top_file(top_file)
    
    # Get charges and neutralize
    # We use MDAnalysis to quickly count atoms to be safe
    u = mda.Universe(gro_file)
    n_atoms = len(u.atoms)
    
    raw_charges = assign_foyer_charges(top_lines, n_atoms)
    clean_charges = neutralize_charges(raw_charges)
    
    # Write Final ITP/TOP
    final_itp = os.path.join(output_dir, f"{resname}.itp")
    write_neutralized_topology(top_lines, final_itp, clean_charges, n_atoms)
    
    # Create a clean master topology file
    final_top = os.path.join(output_dir, "topol.top")
    with open(final_top, 'w') as f:
        f.write(f"; OPLS-AA Topology for {job_name}\n")
        f.write(f'#include "{os.path.basename(ff_xml)}"\n') # Reference the standard OPLS if needed, or include atoms directly
        # For standalone, Foyer puts everything in one, but we split it.
        # Simplified for tutorial: Just point to the ITP
        f.write(f'#include "{os.path.basename(final_itp)}"\n\n')
        f.write("[ system ]\n")
        f.write(f"{job_name}\n\n")
        f.write("[ molecules ]\n")
        f.write(f"Other 1\n") # Foyer often names the residue 'Other' or 'RES'

    print(f"✅ OPLS System Ready!")
    print(f"   -> Topology: {final_itp}")
    print(f"   -> Structure: {gro_file}")
    
    return gro_file, final_itp, final_top
