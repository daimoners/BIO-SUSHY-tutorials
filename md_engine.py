# md_engine.py
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import os
import time

# Try importing ASE (Atomic Simulation Environment)
try:
    from ase import Atoms, units
    from ase.md.langevin import Langevin
    from ase.io import write
    from ase.calculators.calculator import Calculator, all_changes
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False

# --- CUSTOM RDKIT CALCULATOR FOR ASE ---
class RDKitCalculator(Calculator):
    """
    Bridge: Uses RDKit's MMFF force field to calculate Energy & Forces for ASE.
    """
    implemented_properties = ['energy', 'forces']

    def __init__(self, rdkit_mol):
        Calculator.__init__(self)
        self.mol = Chem.Mol(rdkit_mol) # Local copy
        
        # Setup Force Field
        mp = AllChem.MMFFGetMoleculeProperties(self.mol)
        if mp is None: 
            raise ValueError("Force field initialization failed.")
        self.ff = AllChem.MMFFGetMoleculeForceField(self.mol, mp)

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        
        # Sync ASE positions -> RDKit
        conf = self.mol.GetConformer()
        positions = self.atoms.get_positions()
        for i, pos in enumerate(positions):
            conf.SetAtomPosition(i, pos)
            
        # Calc Energy
        if 'energy' in properties:
            e = self.ff.CalcEnergy()
            self.results['energy'] = e * units.kcal / units.mol
            
        # Calc Forces
        if 'forces' in properties:
            grads = self.ff.CalcGrad()
            # Force = -Gradient
            self.results['forces'] = -np.array(grads).reshape((-1, 3)) * units.kcal / units.mol

# --- MAIN SIMULATION FUNCTION ---
def run_nvt_simulation(pdb_file, temp_k=300, steps=1000, dt_fs=1.0, friction=0.02):
    """
    Runs NVT (Constant Number, Volume, Temperature) Dynamics.
    """
    if not ASE_AVAILABLE:
        print("❌ Error: ASE is not installed.")
        return None
    
    print(f"🚀 MD ENGINE: Initializing Simulation ({steps} steps @ {temp_k}K)...")
    
    # 1. Load System
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if not mol: raise ValueError(f"Could not load {pdb_file}")
    
    # 2. Setup Calculator (The Physics)
    calc = RDKitCalculator(mol)
    
    # 3. Convert to ASE Atoms
    pos = mol.GetConformer().GetPositions()
    sym = [a.GetSymbol() for a in mol.GetAtoms()]
    atoms = Atoms(symbols=sym, positions=pos)
    atoms.calc = calc
    
    # 4. Setup Integrator (Langevin Thermostat)
    dyn = Langevin(atoms, timestep=dt_fs*units.fs, temperature_K=temp_k, friction=friction)
    
    # 5. Trajectory Recorder
    output_name = pdb_file.replace(".pdb", "_traj.pdb")
    if os.path.exists(output_name): os.remove(output_name)
    
    print(f"   -> Output file: {output_name}")
    
    def write_frame():
        # Append frame to PDB
        with open(output_name, 'a') as f:
            write(f, atoms, format='pdb')
            
    dyn.attach(write_frame, interval=10) # Save every 10 steps
    
    # 6. RUN
    start_time = time.time()
    dyn.run(steps)
    end_time = time.time()
    
    print(f"✅ SIMULATION COMPLETE in {end_time - start_time:.2f} seconds.")
    return output_name
