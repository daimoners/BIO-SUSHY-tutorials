# Import modules
from rdkit import Chem
from rdkit.Chem import AllChem
import MDAnalysis as mda
import numpy as np
import random
import os
import glob

# --- HELPER FUNCTIONS ---

def estimate_number_of_monomers(input_params: dict, target_num_atoms=220) -> int:
    mol_A = Chem.MolFromSmiles(input_params["smile_rep_unit_A"])
    if mol_A is None: return 0
    num_atoms_A = Chem.rdmolops.AddHs(mol_A).GetNumAtoms() - 2
    
    if input_params["polymer_type"] == "homopolymer":
        return int(target_num_atoms / num_atoms_A)
    else:
        mol_B = Chem.MolFromSmiles(input_params["smile_rep_unit_B"])
        if mol_B is None: return 0
        num_atoms_B = Chem.rdmolops.AddHs(mol_B).GetNumAtoms() - 2
        ratio = input_params.get("stoichiometric_ratio", 0.5)
        avg_atoms = ratio * num_atoms_A + (1 - ratio) * num_atoms_B
        return int(target_num_atoms / avg_atoms)

# --- MAIN CLASS ---

class PolymerBuilder:
    def __init__(self, smiles_A: str, n_chains: int, polymer_type: str, smile_B="", stoichiometric_ratio=1):
        self.n_chains = n_chains
        self.num_A = int(n_chains * stoichiometric_ratio)
        self.num_B = self.n_chains - self.num_A
        self.res_seq = self.define_res_seq(polymer_type)

        # Setup Heads/Tails
        if self.res_seq[0] == "A": head_smile = smiles_A.replace("[*]", "", 1)
        else: head_smile = smile_B.replace("[*]", "", 1)

        if self.res_seq[-1] == "A": tail_smile = smiles_A[::-1].replace("]*[", "", 1)[::-1]
        else: tail_smile = smile_B[::-1].replace("]*[", "", 1)[::-1]

        self.head = Chem.rdmolops.AddHs(Chem.MolFromSmiles(head_smile))
        self.center_A = Chem.rdmolops.AddHs(Chem.MolFromSmiles(smiles_A))
        self.center_B = Chem.rdmolops.AddHs(Chem.MolFromSmiles(smile_B)) if smile_B else None
        self.tail = Chem.rdmolops.AddHs(Chem.MolFromSmiles(tail_smile))
        
        self.char_2_smile = { "A" : self.center_A, "B" : self.center_B }
        self.polymer_chain = self.head.__copy__()

    def define_res_seq(self, polymer_type:str) -> str:
        if polymer_type == "homopolymer": return "A" * self.n_chains
        elif polymer_type == "copolymer-alternate":
            res_seq = "AB" * int(self.n_chains/2)
            if self.n_chains % 2 != 0: res_seq += "A"
            return res_seq
        elif polymer_type == "copolymer-block": return "A" * self.num_A + "B" * self.num_B
        elif polymer_type == "copolymer-random":
            seq = ["A"] * self.num_A + ["B"] * self.num_B
            random.shuffle(seq)
            return "".join(seq)
        return "A" * self.n_chains

    def add_monomer(self, monomer: Chem.Mol) -> None:
        self.polymer_chain = Chem.rdmolops.CombineMols(self.polymer_chain, monomer)
        dummy_atoms = [atom.GetIdx() for atom in self.polymer_chain.GetAtoms() if atom.GetSymbol() == "*"]
        if len(dummy_atoms) < 2: return 
        
        d1, d2 = dummy_atoms[0], dummy_atoms[1]
        l1 = self.polymer_chain.GetAtomWithIdx(d1).GetNeighbors()[0]
        l2 = self.polymer_chain.GetAtomWithIdx(d2).GetNeighbors()[0]

        ed = Chem.EditableMol(self.polymer_chain)
        ed.RemoveBond(d1, l1.GetIdx()); ed.RemoveBond(d2, l2.GetIdx())
        ed.AddBond(l1.GetIdx(), l2.GetIdx(), Chem.BondType.SINGLE)
        ed.RemoveAtom(max(d1, d2)); ed.RemoveAtom(min(d1, d2))
        self.polymer_chain = ed.GetMol()

    def build_raw_3d(self) -> None:
        print(f"   -> Assembling {self.n_chains} monomers...")
        for char in self.res_seq[1:-1:]:
            self.add_monomer(self.char_2_smile[char])
        self.add_monomer(self.tail)
        Chem.SanitizeMol(self.polymer_chain)
        print("   -> Generating initial 3D coordinates (Embedding)...")
        Chem.AllChem.EmbedMolecule(self.polymer_chain, useRandomCoords=True, randomSeed=42)

    def save_pdb(self, filename: str, resname="POL"):
        mol_tmp = "temp.pdb"
        Chem.MolToPDBFile(self.polymer_chain, mol_tmp)
        u = mda.Universe(mol_tmp)
        u.add_TopologyAttr("resname", [resname])
        
        pos = u.atoms.positions
        u.atoms.positions -= (np.min(pos, axis=0) - 10.0)
        L = np.max(u.atoms.positions) + 10.0
        u.dimensions = [L, L, L, 90, 90, 90]
        
        u.atoms.write(filename)
        os.remove(mol_tmp)
        print(f"   -> Saved: {filename}")

# --- COLAB INTERFACE FUNCTIONS ---

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
        print(f"❌ Error: Invalid SMILES for {label}")
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

def build_polymer(smiles_A, polymer_type, target_size, smiles_B="", ratio=0.5, name="polymer"):
    """ Step 2: Build Raw Structure """
    config = { "smile_rep_unit_A": smiles_A, "smile_rep_unit_B": smiles_B,
               "polymer_type": polymer_type, "stoichiometric_ratio": ratio }
    
    n_chains = estimate_number_of_monomers(config, target_num_atoms=target_size)
    print(f"1. Target: {target_size} atoms (~{n_chains} monomers)")
    
    builder = PolymerBuilder(smiles_A, n_chains, polymer_type, smiles_B, ratio)
    builder.build_raw_3d() 
    
    # Path handling: 'name' might be a full path (e.g., "PTFE/PTFE")
    # PolymerBuilder.save_pdb works with full paths, so we just append suffix
    output_pdb = f"{name}_raw.pdb"
    builder.save_pdb(output_pdb)
    return output_pdb

def minimize_polymer(input_pdb, name="polymer"):
    """ Step 3: Minimize Existing PDB """
    output_pdb = f"{name}_relaxed.pdb"
    print(f"1. Loading {input_pdb}...")
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None: raise ValueError("Could not read PDB file.")
        
    print("2. Minimizing Energy (Relaxing)...")
    try:
        Chem.AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    except:
        Chem.AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
        
    Chem.MolToPDBFile(mol, "temp_min.pdb")
    u = mda.Universe("temp_min.pdb")
    u.add_TopologyAttr("resname", ["POL"])
    
    pos = u.atoms.positions
    u.atoms.positions -= (np.min(pos, axis=0) - 10.0)
    L = np.max(u.atoms.positions) + 10.0
    u.dimensions = [L, L, L, 90, 90, 90]
    
    u.atoms.write(output_pdb)
    os.remove("temp_min.pdb")
    print(f"   -> Saved: {output_pdb}")
    return output_pdb

def get_relaxed_files():
    # Return all relaxed PDBs in current AND subdirectories
    return sorted(glob.glob("**/*_relaxed.pdb", recursive=True))
