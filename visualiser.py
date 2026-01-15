import os
import py3Dmol

# ==========================================
#        STYLE SETTINGS (Unified)
# ==========================================
def apply_custom_style(view):
    """
    Applies the lab's consistent styling to any py3Dmol view.
    Used for both static monomers and dynamic trajectories.
    """
    
    # 1. Base Style: Clear everything and set Gray Sticks
    # This ensures we don't have leftover 'lines' or 'cartoons'
    view.setStyle({'stick': {'color': '#999999', 'radius': 0.15}})

    # 2. Atoms as Colored Spheres (Overlay)
    # We use 'addStyle' to put spheres ON TOP of the sticks
    s_scale = 0.25
    
    # Standard Organic Elements
    view.addStyle({'elem': 'C'}, {'sphere': {'color': '#32CD32', 'scale': s_scale}}) # LimeGreen
    view.addStyle({'elem': 'H'}, {'sphere': {'color': 'white',   'scale': s_scale}})
    view.addStyle({'elem': 'O'}, {'sphere': {'color': 'red',     'scale': s_scale}})
    view.addStyle({'elem': 'N'}, {'sphere': {'color': 'blue',    'scale': s_scale}})
    view.addStyle({'elem': 'S'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    view.addStyle({'elem': 'P'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    
    # Halogens
    view.addStyle({'elem': 'F'},  {'sphere': {'color': 'cyan',    'scale': s_scale}})
    view.addStyle({'elem': 'Cl'}, {'sphere': {'color': '#00FF00', 'scale': s_scale}})
    view.addStyle({'elem': 'Br'}, {'sphere': {'color': 'brown',   'scale': s_scale}})
    view.addStyle({'elem': 'I'},  {'sphere': {'color': 'purple',  'scale': s_scale}})
    
    # 3. Special/Dummy Atoms (Linkers)
    view.addStyle({'elem': 'X'}, {'sphere': {'color': 'pink',    'scale': 0.4}})
    view.addStyle({'elem': 'R'}, {'sphere': {'color': 'pink',    'scale': 0.4}})

# ==========================================
#        VISUALIZATION FUNCTIONS
# ==========================================

def show_molecule(file_path, width=800, height=400):
    """
    Visualizes a static PDB or MOL file.
    """
    if not file_path or not os.path.exists(file_path):
        print(f"⚠️ File not found: {file_path}")
        return

    view = py3Dmol.view(width=width, height=height)

    # Determine format
    ext = file_path.split(".")[-1].lower()
    with open(file_path, 'r') as f:
        data = f.read()

    if 'mol' in ext:
        view.addModel(data, 'mol')
    else:
        # Default to PDB for everything else (gro, pdb, ent)
        view.addModel(data, 'pdb')

    apply_custom_style(view)
    view.zoomTo()
    view.show()

def show_trajectory(pdb_content, width=800, height=400):
    """
    Visualizes an animated trajectory.
    Args:
        pdb_content (str): The multi-frame PDB data string.
    """
    view = py3Dmol.view(width=width, height=height)
    
    # CRITICAL FIX: explicitly pass 'pdb' as the format.
    # Without this, py3Dmol might not parse the Element column correctly,
    # causing the 'elem' selectors in the style to fail.
    view.addModelsAsFrames(pdb_content, 'pdb')
    
    # Apply the exact same style as the monomer
    apply_custom_style(view)
    
    # Animate
    # loop: forward, backAndForth, or none
    view.animate({'loop': 'forward', 'reps': 50, 'step': 1, 'interval': 60})
    
    view.zoomTo()
    view.show()
