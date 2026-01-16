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
    view.setStyle({'stick': {'color': '#999999', 'radius': 0.15}})

    # 2. Atoms as Colored Spheres (Overlay)
    s_scale = 0.25
    
    # Standard Organic Elements
    view.addStyle({'elem': 'C'}, {'sphere': {'color': '#32CD32', 'scale': s_scale}}) # LimeGreen
    view.addStyle({'elem': 'H'}, {'sphere': {'color': 'white',   'scale': s_scale}})
    view.addStyle({'elem': 'O'}, {'sphere': {'color': 'red',     'scale': s_scale}})
    view.addStyle({'elem': 'N'}, {'sphere': {'color': 'blue',    'scale': s_scale}})
    view.addStyle({'elem': 'S'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    view.addStyle({'elem': 'P'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    
    # Halogens
    view.addStyle({'elem': 'F'},  {'sphere': {'color': 'yellow',    'scale': s_scale}})
    view.addStyle({'elem': 'Cl'}, {'sphere': {'color': '#00FF00', 'scale': s_scale}})
    view.addStyle({'elem': 'Br'}, {'sphere': {'color': 'brown',   'scale': s_scale}})
    view.addStyle({'elem': 'I'},  {'sphere': {'color': 'purple',  'scale': s_scale}})
    
    # 3. Special/Dummy Atoms (Linkers)
    view.addStyle({'elem': 'X'}, {'sphere': {'color': 'pink',    'scale': 0.4}})
    view.addStyle({'elem': 'R'}, {'sphere': {'color': 'pink',    'scale': 0.4}})

# ==========================================
#        VISUALIZATION FUNCTIONS
# ==========================================

def show_molecule(file_path_or_string, width=800, height=400):
    """
    Visualizes a static PDB or MOL file.
    Returns the view object (does not show it immediately).
    """
    view = py3Dmol.view(width=width, height=height)

    # Check if input is a file path or raw string content
    if os.path.exists(file_path_or_string):
        ext = file_path_or_string.split(".")[-1].lower()
        with open(file_path_or_string, 'r') as f:
            data = f.read()
    else:
        # Assume it is raw content if file not found
        data = file_path_or_string
        ext = 'pdb' # Default

    if 'mol' in ext:
        view.addModel(data, 'mol')
    else:
        view.addModel(data, 'pdb')

    apply_custom_style(view)
    view.zoomTo()
    
    return view  # <--- FIXED: Return object instead of showing

def show_trajectory(pdb_content, width=800, height=400):
    """
    Visualizes an animated trajectory.
    Args:
        pdb_content (str): The multi-frame PDB data string.
    Returns:
        view (py3Dmol.view): The view object.
    """
    view = py3Dmol.view(width=width, height=height)
    
    # CRITICAL FIX: explicitly pass 'pdb' as the format.
    view.addModelsAsFrames(pdb_content, 'pdb')
    
    # Apply the exact same style as the monomer
    apply_custom_style(view)
    
    # Animate
    view.animate({'loop': 'forward', 'reps': 50, 'step': 1, 'interval': 60})
    
    view.zoomTo()
    
    return view # <--- FIXED: Return object instead of showing
