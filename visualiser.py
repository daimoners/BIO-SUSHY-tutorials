import os
import py3Dmol

# ==========================================
#        STYLE SETTINGS (Unified)
# ==========================================
def apply_custom_style(view):
    """Applies the lab's consistent styling to any py3Dmol view."""
    
    # 1. Bonds as Gray Sticks
    view.setStyle({'stick': {'color': '#999999', 'radius': 0.15}})

    # 2. Atoms as Colored Spheres
    s_scale = 0.25
    view.addStyle({'elem': 'C'}, {'sphere': {'color': '#32CD32', 'scale': s_scale}}) # LimeGreen
    view.addStyle({'elem': 'H'}, {'sphere': {'color': 'white',   'scale': s_scale}})
    view.addStyle({'elem': 'O'}, {'sphere': {'color': 'red',     'scale': s_scale}})
    view.addStyle({'elem': 'F'}, {'sphere': {'color': 'yellow',  'scale': s_scale}})
    view.addStyle({'elem': 'N'}, {'sphere': {'color': 'blue',    'scale': s_scale}})
    view.addStyle({'elem': 'S'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    view.addStyle({'elem': 'Cl'}, {'sphere': {'color': '#00FF00', 'scale': s_scale}})
    
    # 3. Special Linkers
    view.addStyle({'elem': 'X'}, {'sphere': {'color': 'pink',    'scale': 0.4}})
    view.addStyle({'elem': 'R'}, {'sphere': {'color': 'pink',    'scale': 0.4}})

# ==========================================
#        VISUALIZATION FUNCTIONS
# ==========================================

def show_molecule(file_path, width=800, height=400):
    """Visualizes a static PDB or MOL file."""
    if not file_path or not os.path.exists(file_path):
        print(f"⚠️ File not found: {file_path}")
        return

    ext = file_path.split(".")[-1]
    view = py3Dmol.view(width=width, height=height)

    with open(file_path, 'r') as f:
        data = f.read()

    if ext == 'mol':
        view.addModel(data, 'mol')
    else:
        view.addModel(data, 'pdb')

    apply_custom_style(view)
    view.zoomTo()
    view.show()

def show_trajectory(pdb_content, width=800, height=400):
    """
    Visualizes an animated trajectory.
    Args:
        pdb_content (str): The multi-frame PDB data as a string.
    """
    view = py3Dmol.view(width=width, height=height)
    
    # Load frames
    view.addModelsAsFrames(pdb_content)
    
    # Apply Style
    apply_custom_style(view)
    
    # Animate
    view.animate({'loop': 'forward', 'reps': 0, 'step': 1, 'interval': 50})
    view.zoomTo()
    view.show()
