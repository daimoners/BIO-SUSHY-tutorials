import os
import py3Dmol

# ==========================================
#        STYLING ENGINES
# ==========================================

def style_monomer(view):
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

def style_polymer(view):
    """
    Clean but Chemical: Thick Backbone + Thin Hydrogens.
    Used for: Polymers and Trajectories.
    """
    view.setStyle({}) # Clear previous styles

    # 1. Backbone (Non-Hydrogens): Thicker sticks, colored by element
    #    'greenCarbon' matches the monomer style.
    view.addStyle({'elem': 'H', 'invert': True}, 
                  {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})

    # 2. Hydrogens: Thin white sticks
    #    We bring them back! But radius is small (0.05) so they don't clutter.
    view.addStyle({'elem': 'H'}, 
                  {'stick': {'color': 'white', 'radius': 0.05}})

# ==========================================
#        VISUALIZATION FUNCTIONS
# ==========================================

def show_molecule(file_path, width=800, height=400, style='monomer'):
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

    # Apply Style
    if style == 'polymer':
        style_polymer(view)
    else:
        style_monomer(view)

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
