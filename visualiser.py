import os
import re
import py3Dmol


# ==========================================
#        STYLE SETTINGS (Unified)
# ==========================================

def apply_custom_style(view):
    """
    Applies the lab's consistent styling to any py3Dmol view.
    Used for both static structures and dynamic trajectories.
    """

    # Base style: gray sticks for all bonds.
    view.setStyle({'stick': {'color': '#999999', 'radius': 0.15}})

    # Atom spheres.
    s_scale = 0.25

    # Standard organic elements.
    view.addStyle({'elem': 'C'},  {'sphere': {'color': '#32CD32', 'scale': s_scale}})
    view.addStyle({'elem': 'H'},  {'sphere': {'color': '#E6E6E6', 'scale': s_scale}})
    view.addStyle({'elem': 'O'},  {'sphere': {'color': 'red',     'scale': s_scale}})
    view.addStyle({'elem': 'N'},  {'sphere': {'color': 'blue',    'scale': s_scale}})
    view.addStyle({'elem': 'S'},  {'sphere': {'color': 'orange',  'scale': s_scale}})
    view.addStyle({'elem': 'P'},  {'sphere': {'color': 'orange',  'scale': s_scale}})

    # Halogens.
    view.addStyle({'elem': 'F'},  {'sphere': {'color': 'yellow',  'scale': s_scale}})
    view.addStyle({'elem': 'Cl'}, {'sphere': {'color': '#00FF00', 'scale': s_scale}})
    view.addStyle({'elem': 'Br'}, {'sphere': {'color': 'brown',   'scale': s_scale}})
    view.addStyle({'elem': 'I'},  {'sphere': {'color': 'purple',  'scale': s_scale}})

    # Dummy/linker atoms if present.
    view.addStyle({'elem': 'X'}, {'sphere': {'color': 'pink', 'scale': 0.4}})
    view.addStyle({'elem': 'R'}, {'sphere': {'color': 'pink', 'scale': 0.4}})


# ==========================================
#        INTERNAL HELPERS
# ==========================================

def _guess_element(atom_name):
    """
    Guess a chemical element from an atom name.

    PDB/GRO files generated automatically can have missing or inconsistent
    element columns. This simple guesser is robust enough for the tutorial's
    organic/polymer systems.
    """
    name = (atom_name or "").strip()
    if not name:
        return "C"

    letters = re.sub(r"[^A-Za-z]", "", name)
    if not letters:
        return "C"

    upper = letters.upper()

    if upper.startswith("CL"):
        return "Cl"
    if upper.startswith("BR"):
        return "Br"

    first = upper[0]
    if first in ["H", "C", "N", "O", "S", "P", "F", "I"]:
        return first

    return "C"


def _normalize_pdb_block(pdb_text):
    """
    Rewrite ATOM/HETATM lines so the element field is populated consistently.
    """
    out_lines = []

    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16].strip()
            element = _guess_element(atom_name)
            padded = line.ljust(78)
            fixed = padded[:76] + f"{element:>2s}" + padded[78:]
            out_lines.append(fixed)
        else:
            out_lines.append(line)

    return "\n".join(out_lines)


def _read_gro_file(gro_file):
    """Read a GRO file and return atom records plus orthorhombic box lengths."""
    with open(gro_file, "r", encoding="utf-8", errors="ignore") as handle:
        lines = handle.read().splitlines()

    if len(lines) < 3:
        raise ValueError(f"{gro_file} does not look like a valid GRO file.")

    try:
        n_atoms = int(lines[1].strip())
    except Exception as exc:
        raise ValueError(f"Could not read atom count from {gro_file}.") from exc

    atom_lines = lines[2:2 + n_atoms]
    if len(atom_lines) != n_atoms:
        raise ValueError(f"{gro_file} ended before all atoms were read.")

    box_index = 2 + n_atoms
    if box_index >= len(lines):
        raise ValueError(f"{gro_file} does not contain a box line.")

    box_parts = lines[box_index].split()
    if len(box_parts) < 3:
        raise ValueError(f"{gro_file} box line must contain at least three values.")

    box_nm = [float(x) for x in box_parts[:3]]

    atoms = []
    for line in atom_lines:
        atoms.append({
            "resnr": int(line[0:5]),
            "resname": line[5:10].strip() or "POL",
            "atomname": line[10:15].strip() or "C",
            "atomnr": int(line[15:20]),
            "x_nm": float(line[20:28]),
            "y_nm": float(line[28:36]),
            "z_nm": float(line[36:44]),
        })

    return atoms, box_nm


def _gro_to_pdb_block(gro_file):
    """Convert GRO to a minimal PDB block for py3Dmol."""
    atoms, box_nm = _read_gro_file(gro_file)

    lx, ly, lz = [v * 10.0 for v in box_nm]  # nm -> Angstrom

    pdb_lines = [
        f"CRYST1{lx:9.3f}{ly:9.3f}{lz:9.3f}"
        f"{90.00:7.2f}{90.00:7.2f}{90.00:7.2f} P 1           1"
    ]

    for serial, atom in enumerate(atoms, start=1):
        atomname = atom["atomname"]
        resname = atom["resname"]
        resnr = atom["resnr"]
        x = atom["x_nm"] * 10.0
        y = atom["y_nm"] * 10.0
        z = atom["z_nm"] * 10.0
        element = _guess_element(atomname)

        pdb_lines.append(
            f"ATOM  {serial:5d} {atomname[:4]:>4s} {resname[:3]:>3s} A{resnr % 10000:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {element:>2s}"
        )

    pdb_lines.append("END")

    return "\n".join(pdb_lines), box_nm


def _add_periodic_box_edges(view, box_nm, color="black"):
    """Draw the edges of an orthorhombic periodic box in py3Dmol."""
    if box_nm is None:
        return

    lx, ly, lz = [float(v) * 10.0 for v in box_nm]  # nm -> Angstrom

    corners = [
        (0.0, 0.0, 0.0),
        (lx, 0.0, 0.0),
        (0.0, ly, 0.0),
        (0.0, 0.0, lz),
        (lx, ly, 0.0),
        (lx, 0.0, lz),
        (0.0, ly, lz),
        (lx, ly, lz),
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
            "start": {"x": x1, "y": y1, "z": z1},
            "end": {"x": x2, "y": y2, "z": z2},
            "color": color,
            "linewidth": 2,
        })

def _parse_cryst1_box_nm(pdb_text):
    """
    Parse the first CRYST1 record found in a PDB block.
    Returns box lengths in nm as (lx, ly, lz), or None if missing.
    Assumes orthorhombic boxes.
    """
    for line in pdb_text.splitlines():
        if line.startswith("CRYST1"):
            try:
                a = float(line[6:15])   # Angstrom
                b = float(line[15:24])  # Angstrom
                c = float(line[24:33])  # Angstrom
                alpha = float(line[33:40])
                beta = float(line[40:47])
                gamma = float(line[47:54])

                # This tutorial assumes orthorhombic boxes.
                if abs(alpha - 90.0) < 1e-3 and abs(beta - 90.0) < 1e-3 and abs(gamma - 90.0) < 1e-3:
                    return (a / 10.0, b / 10.0, c / 10.0)  # Å -> nm
            except Exception:
                return None

    return None
    
# ==========================================
#        VISUALIZATION FUNCTIONS
# ==========================================

def show_structure(file_path_or_string, width=800, height=400, show_box="auto"):
    """
    Visualize a static molecular structure using one unified style.

    Supports:
    - MOL files
    - PDB files
    - GRO files
    - raw PDB/MOL strings

    Parameters
    ----------
    file_path_or_string : str
        Path to a .mol, .pdb or .gro file, or raw molecular text.
    width, height : int
        Viewer dimensions.
    show_box : bool or "auto"
        If True, draw the periodic box when available.
        If False, never draw it.
        If "auto", draw it automatically for GRO files.

    Returns
    -------
    py3Dmol.view
    """
    view = py3Dmol.view(width=width, height=height)

    box_nm = None
    ext = "pdb"

    if os.path.exists(file_path_or_string):
        ext = os.path.splitext(file_path_or_string)[1].lower().replace(".", "")

        if ext == "gro":
            data, box_nm = _gro_to_pdb_block(file_path_or_string)
            fmt = "pdb"

        elif ext == "mol":
            with open(file_path_or_string, "r") as handle:
                data = handle.read()
            fmt = "mol"

        else:
            with open(file_path_or_string, "r") as handle:
                data = handle.read()
            data = _normalize_pdb_block(data)
            box_nm = _parse_cryst1_box_nm(data)
            fmt = "pdb"

    else:
        # Raw molecular string. Treat as PDB by default.
        data = _normalize_pdb_block(file_path_or_string)
        box_nm = _parse_cryst1_box_nm(data)
        fmt = "pdb"

    view.addModel(data, fmt)
    apply_custom_style(view)

    if show_box is True or (show_box == "auto" and ext == "gro"):
        if box_nm is not None:
            _add_periodic_box_edges(view, box_nm)

    view.zoomTo()
    return view


def show_trajectory(pdb_content, width=800, height=400, show_box=False):
    """
    Visualize an animated trajectory from a multi-frame PDB string.

    Parameters
    ----------
    pdb_content : str
        Multi-frame PDB content.
    width, height : int
        Viewer dimensions.
    show_box : bool
        If True, try to show the periodic box from CRYST1 records.

    Returns
    -------
    py3Dmol.view
    """
    view = py3Dmol.view(width=width, height=height)

    pdb_content = _normalize_pdb_block(pdb_content)

    view.addModelsAsFrames(pdb_content, "pdb")
    apply_custom_style(view)

    if show_box:
        # Preferred route: use the unit cell embedded in PDB CRYST1 records.
        try:
            view.addUnitCell()
        except Exception:
            # Fallback: draw a static box from the first CRYST1 record.
            box_nm = _parse_cryst1_box_nm(pdb_content)
            if box_nm is not None:
                _add_periodic_box_edges(view, box_nm)

    view.animate({
        "loop": "forward",
        "reps": 50,
        "step": 1,
        "interval": 60,
    })

    view.zoomTo()
    return view
