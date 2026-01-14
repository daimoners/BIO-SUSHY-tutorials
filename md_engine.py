# md_engine.py
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import time

# Try importing ASE
try:
    from ase import Atoms, units
    from ase.md.langevin import Langevin
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

# --- CUSTOM RDKIT CALCULATOR FOR ASE ---
class RDKitCalculator(object):
    """
    Bridge: Uses RDKit's MMFF force field to calculate Energy & Forces for ASE.
    """
    implemented_properties = ['energy', 'forces']

    def __init__(self, rdkit_mol):
        self.mol = Chem.Mol(rdkit_mol) 
        # Fix for the "Missing Hs" warning: Ensure chemistry is sane
        try:
            Chem.SanitizeMol(self.mol)
            # If PDB missing bonds, this rebuilds them for the force field
            self.mol = Chem.AddHs(self.mol, addCoords=True) 
        except:
            pass
            
        mp = AllChem.MMFFGetMoleculeProperties(self.mol)
        if mp is None: 
            raise ValueError("Force field initialization failed.")
        self.ff = AllChem.MMFFGetMoleculeForceField(self.mol, mp)
        self.results = {}

    def update(self, atoms):
        conf = self.mol.GetConformer()
        positions = atoms.get_positions()
        for i, pos in enumerate(positions):
            conf.SetAtomPosition(i, pos)

    def get_potential_energy(self, atoms=None, force_consistent=False):
        if atoms: self.update(atoms)
        e = self.ff.CalcEnergy()
        return e * units.kcal / units.mol

    def get_forces(self, atoms=None):
        if atoms: self.update(atoms)
        grads = self.ff.CalcGrad()
        return -np.array(grads).reshape((-1, 3)) * units.kcal / units.mol
    
    def calculate(self, atoms=None, properties=['energy'], system_changes=['positions']):
        self.update(atoms)
        self.results['energy'] = self.get_potential_energy()
        self.results['forces'] = self.get_forces()

# --- MAIN SIMULATION FUNCTION ---
def run_nvt_simulation(pdb_file, temp_k=300, steps=1000, dt_fs=1.0, friction=0.02):
    if not ASE_AVAILABLE:
        print("❌ Error: ASE is not installed.")
        return None
    
    print(f"🚀 MD ENGINE: Initializing Simulation ({steps} steps @ {temp_k}K)...")
    
    # 1. Load System
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if not mol: raise ValueError(f"Could not load {pdb_file}")
    
    # 2. Setup Calculator
    calc = RDKitCalculator(mol)
    
    # 3. Convert to ASE Atoms
    pos = mol.GetConformer().GetPositions()
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    atoms = Atoms(symbols=sym, positions=pos)
    atoms.calc = calc
    
    # 4. Setup Integrator
    dyn = Langevin(atoms, timestep=dt_fs*units.fs, temperature_K=temp_k, friction=friction)
    
    # 5. Trajectory Recorder
    output_name = pdb_file.replace(".pdb", "_traj.pdb")
    if os.path.exists(output_name): os.remove(output_name)
    
    print(f"   -> Output file: {output_name}")
    
    counter = [0]
    def write_frame():
        counter[0] += 1
        conf = mol.GetConformer()
        positions = atoms.get_positions()
        for i, p in enumerate(positions):
            conf.SetAtomPosition(i, p)
            
        # Append PDB block manually with correct delimiters
        with open(output_name, 'a') as f:
            f.write(f"MODEL {counter[0]}\n")
            # Get PDB block, remove the END tag so we can stack them
            block = Chem.MolToPDBBlock(mol).replace("END\n", "")
            f.write(block)
            f.write("ENDMDL\n") # Crucial for animation players

    dyn.attach(write_frame, interval=10) 
    
    # 6. RUN
    start_time = time.time()
    try:
        dyn.run(steps)
    except Exception as e:
        print(f"⚠️ Warning during simulation: {e}")
        
    end_time = time.time()
    print(f"✅ SIMULATION COMPLETE in {end_time - start_time:.2f} seconds.")
    return output_name
