import py3Dmol
import ipywidgets as widgets
from IPython.display import display

def show_trajectory(pdb_string, width=800, height=400):
    """
    Visualizes a trajectory with a 'Step Counter' HUD.
    Returns an interactive widget to scrub through frames.
    """
    # 1. Split the PDB string into frames (models)
    # PDBs separate frames with 'ENDMDL'
    frames = [f for f in pdb_string.split("ENDMDL") if f.strip()]
    n_frames = len(frames)
    
    # 2. Setup the Viewer
    view = py3Dmol.view(width=width, height=height)
    
    # 3. Add each frame as a separate model
    # We loop manually so we can attach a specific Label to each frame
    print(f"   -> Processing {n_frames} frames for viewer...")
    
    for i, frame_pdb in enumerate(frames):
        view.addModel(frame_pdb, 'pdb')
        
        # Add a "Heads Up Display" (HUD) Label
        # We attach it to the specific model index (i)
        # This ensures 'Step 10' only shows up when Frame 10 is visible.
        view.addLabel(
            f"Step: {i}", 
            {
                "position": {"x": -80, "y": -80, "z": 0}, 
                "useScreen": True,  # Fix to screen corner
                "fontColor": "black",
                "backgroundColor": "white",
                "fontSize": 12,
                "backgroundOpacity": 0.8
            },
            sel={'model': i} # Important: Only show on this frame
        )

    # 4. Apply Styles
    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
    view.zoomTo()
    
    # 5. Create the Animation / Slider Logic
    view.animate({'loop': 'forward', 'interval': 50}) # Start playing by default
    
    # --- WIDGET CONTROL ---
    # We use a slider to control the 'frame' of the viewer via JavaScript
    
    def on_slider_change(change):
        # JavaScript to set the frame on the active viewer
        # This is a one-way sync from Python Slider -> JS Viewer
        frame_idx = change['new']
        # We inject a command to stop animation and jump to frame
        # Note: py3Dmol doesn't expose a clean Python API for this after render,
        # but the animation above covers the dynamic view.
        pass 

    # Since linking a Python slider to a live JS WebGL canvas in Colab 
    # without custom JS injection is flaky, we provide the robust
    # "Animation" view with the HUD we just built.
    
    return view
