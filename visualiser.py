import py3Dmol

def show_trajectory(pdb_string, width=800, height=400):
    """
    Visualizes a trajectory as an animation with a Step Counter HUD.
    """
    view = py3Dmol.view(width=width, height=height)
    
    # 1. Load PDB as Animation Frames
    # This is the key fix: 'addModelsAsFrames' treats the file as a movie
    view.addModelsAsFrames(pdb_string, 'pdb')
    
    # 2. Count frames to add the Step Counter
    # (MDTraj PDBs use 'MODEL' keywords to separate frames)
    count = pdb_string.count("MODEL")
    
    # 3. Add Labels specific to each frame
    # The selector {'model': i} attaches the label ONLY to frame i.
    # When frame i is hidden by the animation, the label hides too.
    print(f"   -> Adding HUD for {count} frames...")
    
    for i in range(count):
        view.addLabel(
            f"Step: {i}", 
            {
                "position": {"x": -80, "y": -80, "z": 0}, 
                "useScreen": True,
                "fontColor": "black",
                "backgroundColor": "white",
                "fontSize": 12,
                "backgroundOpacity": 0.8
            },
            sel={'model': i} # <--- Attaches label to this specific frame
        )

    # 4. Style and Zoom
    view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
    view.zoomTo()
    
    # 5. Play Animation
    view.animate({'loop': 'forward', 'interval': 50})
    
    return view
