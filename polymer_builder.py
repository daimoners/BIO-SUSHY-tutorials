# Import modules
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import glob
import random

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
        
        # Prepare 3D Monomers (Building Blocks)
        # We still embed the *single* monomer to get a valid starting block
        self.block_A = self._prepare_3d_block(smiles_A)
        self.block_B = self._prepare_3d_block(smile_B) if smile_B else None

    def _prepare_3d_block(self, smiles):
        """ Generates a 3D conformer for a single monomer """
        if not smiles: return None
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol) # Optimizing the SINGLE monomer is safe and good
        return mol

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

    def build_linear_chain(self):
        """ 
        Assembles the polymer linearly. 
        NO MINIMIZATION is performed on the full chain to avoid knots.
        """
        print(f"   -> Assembling {self.n_chains} monomers linearly...")
        
        # Start with the first block
        current_type = self.res_seq[0]
        # Copy the monomer
        full_mol = Chem.RWMol(self.block_A if current_type == "A" else self.block_B)
        
        # Shift Vector: 
        # 2.2 Angstroms is a safe "Stretched Bond".
        # It is > 1.54 (C-C bond) to prevent ANY overlap.
        # It is < 3.0 to ensure OpenMM/Gromacs generation doesn't get confused.
        shift_vector = np.array([2.2, 0.0, 0.0]) 
        current_pos_offset = np.array([0.0, 0.0, 0.0])
        
        # Track bonds to create later
        bonds_to_make = [] 
        
        # Identify the first tail (dummy atom)
        last_tail_idx = -1
        for atom in full_mol.GetAtoms():
            if atom.GetSymbol() == "*":
                last_tail_idx = atom.GetIdx() 
                # Pick the first one we find as the tail for now
        
        # Iteratively add monomers
        for i, char in enumerate(self.res_seq[1:]):
            block = self.block_A if char == "A" else self.block_B
            new_monomer = Chem.Mol(block)
            
            # Shift coordinates
            conf = new_monomer.GetConformer(0)
            current_pos_offset += shift_vector
            
            for k in range(new_monomer.GetNumAtoms()):
                pos = conf.GetAtomPosition(k)
                new_pos = pos + current_pos_offset
                conf.SetAtomPosition(k, new_pos)
            
            # Combine
            index_offset = full_mol.GetNumAtoms()
            
            # Find Head of new monomer
            new_head_local_idx = -1
            for atom in new_monomer.GetAtoms():
                if atom.GetSymbol() == "*":
                    new_head_local_idx = atom.GetIdx()
                    break
            
            full_mol = Chem.RWMol(Chem.CombineMols(full_mol, new_monomer))
            
            # Record Bond (Previous Tail -> New Head)
            if last_tail_idx != -1 and new_head_local_idx != -1:
                real_new_head = new_head_local_idx + index_offset
                bonds_to_make.append((last_tail_idx, real_new_head))
                
                # Find the NEW tail for the next iteration
                for idx in range(index_offset, full_mol.GetNumAtoms()):
                    atom = full_mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() == "*" and idx != real_new_head:
                        last_tail_idx = idx
                        break
        
        # Stitch bonds
        print("   -> Stitching bonds...")
        for a1, a2 in bonds_to_make:
            n1 = full_mol.GetAtomWithIdx(a1).GetNeighbors()[0]
            n2 = full_mol.GetAtomWithIdx(a2).GetNeighbors()[0]
            full_mol.AddBond(n1.GetIdx(), n2.GetIdx(), Chem.BondType.SINGLE)
        
        # Remove Dummy Atoms
        atoms_to_remove = [atom.GetIdx() for atom in full_mol.GetAtoms() if atom.GetSymbol() == "*"]
        atoms_to_remove.sort(reverse=True)
        for idx in atoms_to_remove:
            full_mol.RemoveAtom(idx)

        # Finalize
        self.polymer_chain = full_mol.GetMol()
        Chem.SanitizeMol(self.polymer_chain)
        print("   -> Linear chain built (Raw).")

    def save_pdb(self, filename: str, resname="POL"):
        for atom in self.polymer_chain.GetAtoms():
            info = Chem.AtomPDBResidueInfo()
            info.SetResidueName(resname)
            info.SetResidueNumber(1)
            info.SetIsHeteroAtom(False)
            atom.SetPDBResidueInfo(info)
            
        Chem.MolToPDBFile(self.polymer_chain, filename)
        print(f"   -> Saved: {filename}")

# --- COLAB INTERFACE FUNCTIONS ---

def inspect_monomer(smiles, label="monomer"):
    """ Generates a 3D MOL file for a monomer SMILES """
    if not smiles: return None
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print(f"❌ Error: Invalid SMILES for {label}")
        return None

    vis_mol = Chem.RWMol(mol)
    vis_mol = Chem.AddHs(vis_mol) 
    wildcard_idxs = [atom.GetIdx() for atom in vis_mol.GetAtoms() if atom.GetSymbol() == "*"]
    for idx in wildcard_idxs:
        vis_mol.GetAtomWithIdx(idx).SetAtomicNum(6) 
        vis_mol.GetAtomWithIdx(idx).SetHybridization(Chem.rdchem.HybridizationType.SP3)

    Chem.SanitizeMol(vis_mol)
    vis_mol = Chem.AddHs(vis_mol)

    try:
        AllChem.EmbedMolecule(vis_mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(vis_mol)
    except: pass

    directory = os.path.dirname(label)
    if directory: os.makedirs(directory, exist_ok=True)
    filename = label if label.endswith(".mol") else f"{label}.mol"
    
    Chem.MolToMolFile(vis_mol, filename)
    return filename

def build_polymer(smiles_A, polymer_type, target_size, smiles_B="", ratio=0.5, name="polymer"):
    """ Step 2: Build Linear Structure """
    config = { "smile_rep_unit_A": smiles_A, "smile_rep_unit_B": smiles_B,
               "polymer_type": polymer_type, "stoichiometric_ratio": ratio }
    
    n_chains = estimate_number_of_monomers(config, target_num_atoms=target_size)
    print(f"1. Target: {target_size} atoms (~{n_chains} monomers)")
    
    builder = PolymerBuilder(smiles_A, n_chains, polymer_type, smiles_B, ratio)
    builder.build_linear_chain()
    
    output_pdb = f"{name}_polymer.pdb"
    builder.save_pdb(output_pdb)
    return output_pdb

def minimize_polymer(input_pdb, name="polymer"):
    """ 
    Step 3: Finalize Structure (NO MINIMIZATION)
    Just re-saves the file to ensure naming consistency and RDKit parsing.
    The actual energy minimization will happen in MD Engine.
    """
    output_pdb = f"{name}_relaxed.pdb"
    print(f"1. Loading {input_pdb}...")
    
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None: raise ValueError("Could not read PDB file.")
        
    print("2. Formatting PDB (Skipping Minimization)...")
    # We deliberately SKIP Chem.AllChem.UFFOptimizeMolecule here.
    
    Chem.MolToPDBFile(mol, output_pdb)
    print(f"   -> Saved: {output_pdb}")
    return output_pdb

def get_relaxed_files():
    return sorted(glob.glob("**/*_relaxed.pdb", recursive=True))
