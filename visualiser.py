import os
import re
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
#        INTERNAL HELPERS FOR GRO/BULK VIEW
# ==========================================

def _guess_element(atom_name):
    """Guess a chemical element from a PDB/GRO atom name."""
    letters = re.sub(r"[^A-Za-z]", "", atom_name or "").upper()
    if not letters:
        return "C"
    if letters.startswith("CL"):
        return "Cl"
    if letters.startswith("BR"):
        return "Br"
    if letters.startswith("NA"):
        return "Na"
    if letters.startswith("MG"):
        return "Mg"
    if letters.startswith("CA") and atom_name.strip().upper().startswith("CA"):
        # In organic polymer files CA is often an atom name for carbon alpha,
        # but calcium is rare here. Use carbon as the safer tutorial default.
        return "C"
    return letters[0].upper()


def _read_gro_for_visualisation(gro_file):
    """Read atom records and orthorhombic box lengths from a GRO file."""
    with open(gro_file, "r") as handle:
        lines = handle.read().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{gro_file} does not look like a valid GRO file.")

    n_atoms = int(lines[1].strip())
    atom_lines = lines[2:2 + n_atoms]
    box_values = [float(x) for x in lines[2 + n_atoms].split()[:3]]
    if len(box_values) < 3:
        raise ValueError("The GRO box line must contain at least three box lengths.")

    atoms = []
    for line in atom_lines:
        # Standard GRO fixed-width fields. Coordinates are in nm.
        resnr = int(line[0:5])
        resname = line[5:10].strip() or "POL"
        atomname = line[10:15].strip() or "C"
        atomnr = int(line[15:20])
        x_nm = float(line[20:28])
        y_nm = float(line[28:36])
        z_nm = float(line[36:44])
        atoms.append((resnr, resname, atomname, atomnr, x_nm, y_nm, z_nm))
    return atoms, box_values


def _gro_to_pdb_block(gro_file, max_atoms=None):
    """
    Convert a GRO structure to a minimal PDB block for py3Dmol.

    GRO coordinates are in nm, while PDB coordinates are in Angstrom.
    """
    atoms, box_nm = _read_gro_for_visualisation(gro_file)
    original_count = len(atoms)
    if max_atoms is not None and original_count > int(max_atoms):
        stride = max(1, original_count // int(max_atoms))
        atoms = atoms[::stride]

    lx, ly, lz = [v * 10.0 for v in box_nm]  # nm -> Angstrom
    pdb_lines = [
        f"CRYST1{lx:9.3f}{ly:9.3f}{lz:9.3f}{90.00:7.2f}{90.00:7.2f}{90.00:7.2f} P 1           1"
    ]

    for serial, (resnr, resname, atomname, _atomnr, x_nm, y_nm, z_nm) in enumerate(atoms, start=1):
        element = _guess_element(atomname)
        x, y, z = x_nm * 10.0, y_nm * 10.0, z_nm * 10.0
        pdb_lines.append(
            f"ATOM  {serial:5d} {atomname[:4]:>4s} {resname[:3]:>3s} A{resnr % 10000:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}"
        )
    pdb_lines.append("END")
    return "\n".join(pdb_lines), box_nm, original_count, len(atoms)


def _add_periodic_box_edges(view, box_nm, color="black"):
    """Draw the edges of an orthorhombic periodic box in a py3Dmol view."""
    lx, ly, lz = [float(v) * 10.0 for v in box_nm]  # nm -> Angstrom
    corners = [
        (0.0, 0.0, 0.0), (lx, 0.0, 0.0), (0.0, ly, 0.0), (0.0, 0.0, lz),
        (lx, ly, 0.0), (lx, 0.0, lz), (0.0, ly, lz), (lx, ly, lz),
    ]
    edges = [
        (0, 1), (0, 2), (0, 3),
        (1, 4), (1, 5),
        (2, 4), (2, 6),
        (3, 5), (3, 6),
        (4, 7), (5, 7), (6, 7),
    ]
    for i, j in edges:
        x1, y1, z1 = corners[i]
        x2, y2, z2 = corners[j]
        view.addLine({
            'start': {'x': x1, 'y': y1, 'z': z1},
            'end': {'x': x2, 'y': y2, 'z': z2},
            'color': color,
            'linewidth': 2,
        })


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
    elif ext == 'gro':
        data, _, _, _ = _gro_to_pdb_block(file_path_or_string)
        view.addModel(data, 'pdb')
    else:
        view.addModel(data, 'pdb')

    apply_custom_style(view)
    view.zoomTo()
    
    return view  # <--- FIXED: Return object instead of showing


def show_bulk_box(gro_file, width=850, height=500, max_atoms=3000):
    """
    Visualize a mini-bulk GRO file together with its periodic simulation box.

    Parameters
    ----------
    gro_file : str
        Path to a GRO file containing the mini-bulk coordinates and box.
    width, height : int
        Viewer dimensions.
    max_atoms : int or None
        Maximum number of atoms shown in the py3Dmol model. Large systems can
        be slow in a browser; if the structure is larger, atoms are shown with
        a regular stride for visualization only.
    """
    pdb_block, box_nm, original_count, shown_count = _gro_to_pdb_block(gro_file, max_atoms=max_atoms)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_block, 'pdb')
    apply_custom_style(view)
    _add_periodic_box_edges(view, box_nm)
    view.zoomTo()

    if shown_count < original_count:
        print(
            f"⚠️ Visualisation sampled {shown_count}/{original_count} atoms for browser performance. "
            "The simulation still uses the full structure."
        )
    print(f"📦 Box shown: {box_nm[0]:.3f} × {box_nm[1]:.3f} × {box_nm[2]:.3f} nm")
    return view


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
