import os
import py3Dmol

# --- 2. DEFINE CUSTOM VISUALIZER (FIXED FOR MONOMERS) ---
def show_molecule(file_path, width=800, height=400):
    """
    Visualizes PDB (Polymers) or MOL (Monomers) files.
    """
    if not file_path or not os.path.exists(file_path):
        print(f"⚠️ File not found: {file_path}")
        return

    ext = file_path.split(".")[-1] # Get extension (pdb or mol)
    view = py3Dmol.view(width=width, height=height)

    with open(file_path, 'r') as f:
        data = f.read()

    # Load data based on format
    if ext == 'mol':
        view.addModel(data, 'mol') # Strict connectivity
    else:
        view.addModel(data, 'pdb') # Standard PDB

    # 1. Set Bonds to Gray Sticks
    view.setStyle({'stick': {'color': '#999999', 'radius': 0.15}})

    # 2. Add Colored Atoms (Spheres)
    s_scale = 0.25
    view.addStyle({'elem': 'C'}, {'sphere': {'color': '#32CD32', 'scale': s_scale}}) # LimeGreen
    view.addStyle({'elem': 'H'}, {'sphere': {'color': 'white',   'scale': s_scale}})
    view.addStyle({'elem': 'O'}, {'sphere': {'color': 'red',     'scale': s_scale}})
    view.addStyle({'elem': 'F'}, {'sphere': {'color': 'yellow',  'scale': s_scale}})
    view.addStyle({'elem': 'N'}, {'sphere': {'color': 'blue',    'scale': s_scale}})
    view.addStyle({'elem': 'S'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    view.addStyle({'elem': 'Cl'}, {'sphere': {'color': '#00FF00', 'scale': s_scale}})
    
    # 3. Handle Special "Linker" Atoms [*] if they appear (often X or R)
    view.addStyle({'elem': 'X'}, {'sphere': {'color': 'pink',    'scale': 0.4}})
    view.addStyle({'elem': 'R'}, {'sphere': {'color': 'pink',    'scale': 0.4}})

    view.zoomTo()
    view.show()
    
