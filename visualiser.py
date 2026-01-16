import os
import py3Dmol

# ==========================================
#        STYLING ENGINES
# ==========================================

def style_polymer(view):
    """
    High Detail: Sticks + Spheres + Hydrogens.
    Used for: Single Monomers (Step 3).
    """
    # 1. Base: Gray Sticks
    view.setStyle({'stick': {'color': '#999999', 'radius': 0.15}})

    # 2. Overlay: Colored Spheres (The "Ball & Stick" look)
    s_scale = 0.25
    view.addStyle({'elem': 'C'}, {'sphere': {'color': '#32CD32', 'scale': s_scale}})
    view.addStyle({'elem': 'H'}, {'sphere': {'color': 'white',   'scale': s_scale}})
    view.addStyle({'elem': 'O'}, {'sphere': {'color': 'red',     'scale': s_scale}})
    view.addStyle({'elem': 'N'}, {'sphere': {'color': 'blue',    'scale': s_scale}})
    view.addStyle({'elem': 'S'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    
    # Halogens
    view.addStyle({'elem': 'F'},  {'sphere': {'color': 'yellow',    'scale': s_scale}})
    view.addStyle({'elem': 'Cl'}, {'sphere': {'color': '#00FF00', 'scale': s_scale}})

# ==========================================
#        VISUALIZATION FUNCTIONS
# ==========================================

def show_molecule(file_path, width=800, height=400, style='polymer'):
    if not file_path or not os.path.exists(file_path):
        print(f"⚠️ File not found: {file_path}")
        return

    view = py3Dmol.view(width=width, height=height)

    ext = file_path.split(".")[-1].lower()
    with open(file_path, 'r') as f:
        data = f.read()

    if 'mol' in ext:
        view.addModel(data, 'mol')
    else:
        view.addModel(data, 'pdb')

    view.zoomTo()
    view.show()

def show_trajectory(pdb_content, width=800, height=400):
    view = py3Dmol.view(width=width, height=height)
    view.addModelsAsFrames(pdb_content, 'pdb')
    
    # Always use the Clean Polymer style for motion
    style_polymer(view)
    
    view.animate({'loop': 'forward', 'reps': 50, 'step': 1, 'interval': 60})
    view.zoomTo()
    view.show()
