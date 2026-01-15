import os
import glob
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

def inspect_monomer(smiles, label="monomer"):
    """
    Generates a 3D MOL file for a monomer SMILES.
    Handles wildcards [*] by temporarily capping them with Carbon 
    so the 3D embedding engine (UFF) doesn't crash or warn.
    """
    if not smiles:
        return None

    # 1. Create Molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print(f"❌ Invalid SMILES: {smiles}")
        return None

    # 2. Prepare for Visualization (Handle Wildcards)
    # We work on a copy so we don't change the actual chemistry logic, just the picture.
    vis_mol = Chem.RWMol(mol)
    vis_mol = Chem.AddHs(vis_mol) # Add H's first

    # Find wildcards [*] and replace them with Carbon (Atomic Num 6)
    # This tricks the physics engine into treating them as normal atoms.
    wildcard_idxs = [atom.GetIdx() for atom in vis_mol.GetAtoms() if atom.GetSymbol() == "*"]
    
    for idx in wildcard_idxs:
        vis_mol.GetAtomWithIdx(idx).SetAtomicNum(6) # Turn * into C
        vis_mol.GetAtomWithIdx(idx).SetHybridization(Chem.rdchem.HybridizationType.SP3)

    # Re-sanitize and Add Hydrogens to the new "Carbons"
    Chem.SanitizeMol(vis_mol)
    vis_mol = Chem.AddHs(vis_mol)

    # 3. Embed 3D Coordinates
    try:
        # Use random coordinates as a seed to ensure generation works
        AllChem.EmbedMolecule(vis_mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(vis_mol)
    except Exception as e:
        print(f"⚠️ Minor 3D generation warning: {e}")

    # 4. Save to File (Path Safe)
    # We split the label into directory and filename to avoid the "viz_/path/" error
    directory = os.path.dirname(label)
    base_name = os.path.basename(label)
    
    if directory:
        os.makedirs(directory, exist_ok=True)
        filename = os.path.join(directory, f"{base_name}.mol")
    else:
        filename = f"viz_{base_name}.mol"

    print(f"   -> Generated 3D view: {filename}")
    Chem.MolToMolFile(vis_mol, filename)
    
    return filename

def build_polymer(smiles_A, polymer_type="homopolymer", target_size=10, smiles_B=None, name="polymer", ratio=0.5):
    """
    Builds a linear polymer chain using the 'Linear' builder from mBuild (implied logic)
    or RDKit string manipulation if mBuild is not available. 
    
    For this tutorial, we use a robust RDKit string stitching method 
    to ensure it works without complex dependencies.
    """
    
    # Clean inputs
    sA = smiles_A.replace("[*]", "").replace("*", "")
    if smiles_B:
        sB = smiles_B.replace("[*]", "").replace("*", "")
    
    # Logic for Sequence
    chain_smiles = ""
    
    if polymer_type == "homopolymer":
        # A-A-A-A...
        chain_smiles = sA * target_size
        
    elif polymer_type == "copolymer-block":
        # A-A-A...B-B-B...
        nA = int(target_size * ratio)
        nB = target_size - nA
        chain_smiles = (sA * nA) + (sB * nB)
        
    elif polymer_type == "copolymer-alternate":
        # A-B-A-B...
        chain_smiles = (sA + sB) * (target_size // 2)
        
    elif polymer_type == "copolymer-random":
        # A-B-A-A-B...
        import random
        block = []
        nA = int(target_size * ratio)
        nB = target_size - nA
        block = [sA]*nA + [sB]*nB
        random.shuffle(block)
        chain_smiles = "".join(block)

    # Create Molecule
    # We add terminal Hydrogens implicitly by parsing SMILES
    mol = Chem.MolFromSmiles(chain_smiles)
    mol = Chem.AddHs(mol)
    
    # Embed 3D (This takes time for large molecules)
    print(f"   -> Embedding {mol.GetNumAtoms()} atoms (this may take a moment)...")
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    
    # Save
    # Handle path vs filename again
    directory = os.path.dirname(name)
    base_name = os.path.basename(name)
    
    if directory:
        os.makedirs(directory, exist_ok=True)
        final_path = os.path.join(directory, f"{base_name}_relaxed.pdb")
    else:
        final_path = f"{base_name}_relaxed.pdb"

    Chem.MolToPDBFile(mol, final_path)
    
    return final_path

def get_relaxed_files():
    """Returns a list of all relaxed PDB files in the current directory."""
    return glob.glob("*_relaxed.pdb")
