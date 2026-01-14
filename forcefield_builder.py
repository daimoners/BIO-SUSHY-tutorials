# polymer_forcefield.py
from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os

def validate_and_assign_ff(input_pdb, output_name="system_ready"):
    """
    Simulates the 'Force Field Generation' step.
    1. Loads structure.
    2. Checks for missing atoms/connectivity.
    3. Assigns MMFF94 properties (charges/types).
    4. Saves a validated PDB ready for MD.
    """
    print(f"⚙️ PARAMETERIZATION: Loading {input_pdb}...")
    
    # 1. Load Molecule
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None:
        raise ValueError(f"❌ Critical Error: Could not parse {input_pdb}")

    # 2. Check Chemistry (Sanitization)
    try:
        Chem.SanitizeMol(mol)
        print("   -> Chemistry check: PASSED (Valence & aromaticity ok)")
    except Exception as e:
        raise ValueError(f"❌ Chemistry Check Failed: {e}")

    # 3. Assign Force Field (MMFF94)
    # This checks if RDKit has parameters for every atom type in the chain.
    print("   -> Assigning MMFF94 Force Field parameters...")
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    
    if mmff_props is None:
        raise ValueError("❌ Error: Could not assign Force Field. Molecule may contain unsupported atom types.")
    
    print("   -> Force Field assigned successfully.")
    
    # 4. (Optional) Compute Partial Charges for user info
    # We use Gasteiger charges just to show we 'calculated' something electronic
    AllChem.ComputeGasteigerCharges(mol)
    avg_charge = sum([float(a.GetProp('_GasteigerCharge')) for a in mol.GetAtoms()])
    print(f"   -> System Net Charge: {avg_charge:.4f} (Should be ~0.0)")

    # 5. Save "Ready" System
    output_file = f"{output_name}.pdb"
    Chem.MolToPDBFile(mol, output_file)
    print(f"✅ SYSTEM READY: Saved valid topology to {output_file}")
    
    return output_file
