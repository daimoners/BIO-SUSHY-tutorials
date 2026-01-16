import os
import py3Dmol

# ==========================================
#        STYLING ENGINES
# ==========================================

def style_monomer(view):
    """
    Detailed 'Ball-and-Stick' style.
    Best for: Single monomers, small molecules.
    Features: Shows Hydrogens, shows double bonds (via gray sticks), highlights atoms with spheres.
    """
    # 1. Base: Gray Sticks (shows bonds clearly)
    view.setStyle({'stick': {'color': '#999999', 'radius': 0.15}})

    # 2. Overlay: Colored Spheres
    s_scale = 0.25
    
    # Standard Organic Elements
    view.addStyle({'elem': 'C'}, {'sphere': {'color': '#32CD32', 'scale': s_scale}}) # LimeGreen
    view.addStyle({'elem': 'H'}, {'sphere': {'color': 'white',   'scale': s_scale}})
    view.addStyle({'elem': 'O'}, {'sphere': {'color': 'red',     'scale': s_scale}})
    view.addStyle({'elem': 'N'}, {'sphere': {'color': 'blue',    'scale': s_scale}})
    view.addStyle({'elem': 'S'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    view.addStyle({'elem': 'P'}, {'sphere': {'color': 'orange',  'scale': s_scale}})
    
    # Halogens & Others
    view.addStyle({'elem': 'F'},  {'sphere': {'color': 'yellow',    'scale': s_scale}})
    view.addStyle({'elem': 'Cl'}, {'sphere': {'color': '#00FF00', 'scale': s_scale}})
    view.addStyle({'elem': 'Br'}, {'sphere': {'color': 'brown',   'scale': s_scale}})
    view.addStyle({'elem': 'I'},  {'sphere': {'color': 'purple',  'scale': s_scale}})

def style_polymer(view):
    """
    Clean 'Stick' style.
    Best for: Long chains, trajectories, dense systems.
    Features: Hides Hydrogens (reduces noise), no spheres, thicker sticks colored by element.
    """
    view.setStyle({}) # Clear previous styles

    # 1. Select ALL atoms EXCEPT Hydrogen (invert=True)
    #    This removes the "fuzz" of H atoms so you can see the backbone coiling.
    # 2. Style: Thicker sticks, colored by element
    #    'greenCarbon' matches your monomer's green spheres.
    view.addStyle({'elem': 'H', 'invert': True}, 
                  {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.25}})

# ==========================================
#        VISUALIZATION FUNCTIONS
# ==========================================

def show_molecule(file_path, width=800, height=400, style='monomer'):
    """
    Visualizes a static structure.
    Args:
        style (str): 'monomer' (detailed) or 'polymer' (clean).
    """
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

    # Apply selected style
    if style == 'polymer':
        style_polymer(view)
    else:
        style_monomer(view)

    view.zoomTo()
    view.show()

def show_trajectory(pdb_content, width=800, height=400):
    """
    Visualizes an animated trajectory.
    Always uses 'polymer' style for performance and clarity.
    """
    view = py3Dmol.view(width=width, height=height)
    
    # Load frames
    view.addModelsAsFrames(pdb_content, 'pdb')
    
    # Apply Clean Polymer Style
    style_polymer(view)
    
    # Animate
    view.animate({'loop': 'forward', 'reps': 50, 'step': 1, 'interval': 60})
    
    view.zoomTo()
    view.show()
