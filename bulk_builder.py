"""
bulk_builder.py

Tutorial bulk builder for BIO-SUSHY.

The user only chooses the number of polymer chains.  The initial cubic box size
is calculated automatically from the chain mass and an internal dilute packing
estimate.  This creates a reasonable starting point for a short NPT relaxation
without presenting the density as a final predicted material property.
"""

import os
import re
import numpy as np
import parmed as pmd


# ==========================================
#        GRO PARSING / WRITING
# ==========================================

def _read_gro_atoms(gro_file):
    """Read atom records from a GRO file. Coordinates are returned in nm."""
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

    atoms = []
    for line in atom_lines:
        # Standard fixed-width GRO columns.
        atoms.append({
            "resnr": int(line[0:5]),
            "resname": line[5:10].strip() or "POL",
            "atomname": line[10:15].strip() or "C",
            "atomnr": int(line[15:20]),
            "x": float(line[20:28]),
            "y": float(line[28:36]),
            "z": float(line[36:44]),
        })

    return atoms


def _write_gro(atoms, box_nm, output_gro, title="Bulk polymer system"):
    """Write a GRO file. Coordinates and box are expected in nm."""
    with open(output_gro, "w", encoding="utf-8") as handle:
        handle.write(f"{title}\n")
        handle.write(f"{len(atoms):5d}\n")

        for i, atom in enumerate(atoms, start=1):
            resnr = int(atom["resnr"]) % 100000
            atomnr = i % 100000
            if atomnr == 0:
                atomnr = 100000

            handle.write(
                f"{resnr:5d}"
                f"{atom['resname'][:5]:<5s}"
                f"{atom['atomname'][:5]:>5s}"
                f"{atomnr:5d}"
                f"{atom['x']:8.3f}"
                f"{atom['y']:8.3f}"
                f"{atom['z']:8.3f}\n"
            )

        handle.write(f"{box_nm[0]:10.5f}{box_nm[1]:10.5f}{box_nm[2]:10.5f}\n")


# ==========================================
#        MASS / BOX ESTIMATION
# ==========================================

def _guess_element(atom_name):
    """Small helper used only as a fallback for mass estimation."""
    letters = re.sub(r"[^A-Za-z]", "", atom_name or "").upper()
    if not letters:
        return "C"
    if letters.startswith("CL"):
        return "Cl"
    if letters.startswith("BR"):
        return "Br"
    return letters[0]


_FALLBACK_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Br": 79.904,
    "I": 126.904,
}


def _calculate_chain_mass_amu(gro_file):
    """
    Calculate chain mass from a GRO file.

    ParmEd is used first.  If masses are missing, a simple element-based
    fallback is used so the tutorial remains robust in Colab.
    """
    try:
        structure = pmd.load_file(gro_file)
        mass = float(sum(atom.mass for atom in structure.atoms))
        if mass > 0:
            return mass
    except Exception:
        pass

    atoms = _read_gro_atoms(gro_file)
    mass = 0.0
    for atom in atoms:
        elem = _guess_element(atom["atomname"])
        mass += _FALLBACK_MASSES.get(elem, 12.011)

    if mass <= 0:
        raise ValueError("Could not determine chain mass from the GRO file.")

    return float(mass)


def _calculate_auto_box_nm(
    chain_mass_amu,
    n_chains,
    chain_extent_nm,
    initial_packing_density_kg_m3=180.0,
    minimum_box_nm=4.5,
    padding_nm=0.8,
):
    """
    Calculate a cubic initial box size.

    The internal density is a dilute packing estimate, not a final material
    density.  The box is also constrained to be large enough for the polymer
    chain extent and for browser-friendly tutorial systems.
    """
    amu_to_kg = 1.66054e-27

    total_mass_kg = float(chain_mass_amu) * float(n_chains) * amu_to_kg
    volume_m3 = total_mass_kg / float(initial_packing_density_kg_m3)
    volume_nm3 = volume_m3 / 1e-27
    density_box_nm = volume_nm3 ** (1.0 / 3.0)

    # Ensure the initial box is not smaller than the current chain extension.
    extent_box_nm = float(chain_extent_nm) + 2.0 * float(padding_nm)

    # Mild additional scaling with the number of chains to reduce initial clashes.
    grid_dim = int(np.ceil(float(n_chains) ** (1.0 / 3.0)))
    grid_box_nm = max(float(chain_extent_nm) * 0.75, 0.8) * grid_dim

    return float(max(density_box_nm + padding_nm, extent_box_nm, grid_box_nm, minimum_box_nm))


# ==========================================
#        GEOMETRY HELPERS
# ==========================================

def _random_rotation_matrix(rng):
    """Generate a random 3D rotation matrix from a random unit quaternion."""
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q

    return np.array([
        [1 - 2*y*y - 2*z*z,     2*x*y - 2*z*w,     2*x*z + 2*y*w],
        [2*x*y + 2*z*w,         1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w],
        [2*x*z - 2*y*w,         2*y*z + 2*x*w,     1 - 2*x*x - 2*y*y],
    ])


def _grid_centres(n_chains, box_nm, rng):
    """Generate roughly uniform placement centres inside a cubic box."""
    grid_dim = int(np.ceil(float(n_chains) ** (1.0 / 3.0)))
    spacing = float(box_nm) / grid_dim

    centres = []
    for ix in range(grid_dim):
        for iy in range(grid_dim):
            for iz in range(grid_dim):
                centre = np.array([
                    (ix + 0.5) * spacing,
                    (iy + 0.5) * spacing,
                    (iz + 0.5) * spacing,
                ])
                jitter = rng.uniform(-0.15 * spacing, 0.15 * spacing, size=3)
                centres.append((centre + jitter) % box_nm)

    rng.shuffle(centres)
    return centres[:n_chains]


# ==========================================
#        TOPOLOGY UPDATE
# ==========================================

def _copy_and_update_topology(single_chain_top, output_top, n_chains):
    """
    Copy a single-chain Gromacs topology and update [ molecules ].

    The first molecule entry is retained and its count is replaced by n_chains.
    Extra molecule lines are skipped for this tutorial workflow.
    """
    with open(single_chain_top, "r", encoding="utf-8", errors="ignore") as handle:
        lines = handle.readlines()

    output_lines = []
    in_molecules = False
    molecule_line_updated = False
    molecule_name = "POL"

    for line in lines:
        stripped = line.strip()

        if stripped.lower() == "[ molecules ]":
            in_molecules = True
            output_lines.append(line)
            continue

        if in_molecules:
            if stripped.startswith("[") and stripped.endswith("]"):
                in_molecules = False
                output_lines.append(line)
                continue

            if stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if len(parts) >= 1 and not molecule_line_updated:
                    molecule_name = parts[0]
                    output_lines.append(f"{molecule_name:<10s} {int(n_chains)}\n")
                    molecule_line_updated = True
                else:
                    # Skip extra molecule lines in the tutorial topology.
                    continue
            else:
                output_lines.append(line)
        else:
            output_lines.append(line)

    if not molecule_line_updated:
        output_lines.append("\n[ molecules ]\n")
        output_lines.append(f"{molecule_name:<10s} {int(n_chains)}\n")

    with open(output_top, "w", encoding="utf-8") as handle:
        handle.writelines(output_lines)


# ==========================================
#        PUBLIC BUILDER
# ==========================================

def build_polymer_bulk(
    single_chain_gro,
    single_chain_top,
    n_chains=6,
    output_dir="Bulk",
    random_seed=42,
    initial_packing_density_kg_m3=180.0,
    minimum_box_nm=4.5,
    padding_nm=0.8,
):
    """
    Build a tutorial mini-bulk system from one parameterized polymer chain.

    User-facing input:
        n_chains

    Internally calculated:
        box size, from chain mass and a dilute initial packing estimate.

    Parameters
    ----------
    single_chain_gro : str
        GRO file of one parameterized polymer chain.
    single_chain_top : str
        TOP file of one parameterized polymer chain.
    n_chains : int
        Number of polymer chains to place in the box.
    output_dir : str
        Output directory.
    random_seed : int
        Random seed for reproducible placement.
    initial_packing_density_kg_m3 : float
        Internal dilute packing-density estimate used only to calculate the
        initial box size. This is not reported as a final material property.
    minimum_box_nm : float
        Lower bound for the cubic box length.
    padding_nm : float
        Extra padding used in the box-size estimate.

    Returns
    -------
    tuple
        bulk_gro, bulk_top, metadata
    """
    n_chains = int(n_chains)
    if n_chains < 1:
        raise ValueError("n_chains must be at least 1.")

    if not os.path.exists(single_chain_gro):
        raise FileNotFoundError(f"Single-chain GRO not found: {single_chain_gro}")
    if not os.path.exists(single_chain_top):
        raise FileNotFoundError(f"Single-chain TOP not found: {single_chain_top}")

    os.makedirs(output_dir, exist_ok=True)

    bulk_gro = os.path.join(output_dir, "bulk.gro")
    bulk_top = os.path.join(output_dir, "topol.top")

    rng = np.random.default_rng(int(random_seed))

    chain_atoms = _read_gro_atoms(single_chain_gro)
    chain_mass_amu = _calculate_chain_mass_amu(single_chain_gro)

    coords = np.array([[a["x"], a["y"], a["z"]] for a in chain_atoms], dtype=float)
    centre = coords.mean(axis=0)
    coords_centered = coords - centre
    chain_extent_nm = float(np.max(np.ptp(coords, axis=0)))

    box_nm = _calculate_auto_box_nm(
        chain_mass_amu=chain_mass_amu,
        n_chains=n_chains,
        chain_extent_nm=chain_extent_nm,
        initial_packing_density_kg_m3=initial_packing_density_kg_m3,
        minimum_box_nm=minimum_box_nm,
        padding_nm=padding_nm,
    )

    print("🧱 Building tutorial bulk")
    print(f"   Chains:              {n_chains}")
    print(f"   Chain mass:          {chain_mass_amu:.2f} amu")
    print(f"   Chain extent:        {chain_extent_nm:.3f} nm")
    print(f"   Auto box size:       {box_nm:.3f} nm")
    print(f"   Packing reference:   {initial_packing_density_kg_m3:.1f} kg/m³")
    print("   Note: this is a dilute starting box, not a final density prediction.")

    centres = _grid_centres(n_chains, box_nm, rng)
    all_atoms = []

    for chain_idx, translation in enumerate(centres):
        rotation = _random_rotation_matrix(rng)
        new_coords = coords_centered @ rotation.T + translation

        for atom, xyz in zip(chain_atoms, new_coords):
            new_atom = dict(atom)
            new_atom["resnr"] = chain_idx + 1
            new_atom["x"] = float(xyz[0] % box_nm)
            new_atom["y"] = float(xyz[1] % box_nm)
            new_atom["z"] = float(xyz[2] % box_nm)
            all_atoms.append(new_atom)

    _write_gro(
        all_atoms,
        np.array([box_nm, box_nm, box_nm], dtype=float),
        bulk_gro,
        title=f"Bulk polymer system: {n_chains} chains",
    )

    _copy_and_update_topology(
        single_chain_top=single_chain_top,
        output_top=bulk_top,
        n_chains=n_chains,
    )

    metadata = {
        "n_chains": int(n_chains),
        "chain_mass_amu": float(chain_mass_amu),
        "chain_extent_nm": float(chain_extent_nm),
        "auto_box_nm": float(box_nm),
        "initial_packing_density_kg_m3": float(initial_packing_density_kg_m3),
        "minimum_box_nm": float(minimum_box_nm),
        "padding_nm": float(padding_nm),
        "random_seed": int(random_seed),
        "note": (
            "Box size is calculated automatically from n_chains and an internal "
            "dilute packing estimate. This is not a final density result."
        ),
    }

    print("✅ Bulk files written:")
    print(f"   GRO: {bulk_gro}")
    print(f"   TOP: {bulk_top}")

    return bulk_gro, bulk_top, metadata
