"""
bulk_builder.py

Small, dependency-light utilities for building a demonstrative polymer mini-bulk
from a single-chain Gromacs coordinate/topology pair.

The goal is tutorial robustness, not production-quality amorphous-cell packing.
For research calculations, replace this simple packer with a validated packing
workflow such as Packmol, Polymatic, mBuild, or an in-house protocol.
"""

from __future__ import annotations

import math
import os
import random
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np

AMU_TO_G_PER_NM3_AT_DENSITY_1 = 1.66053906660e-3


@dataclass
class GroAtom:
    resnr: int
    resname: str
    atomname: str
    atomnr: int
    xyz_nm: np.ndarray
    tail: str = ""


def _parse_gro_atom_line(line: str) -> GroAtom:
    """Parse a standard .gro atom line using fixed-width fields."""
    try:
        resnr = int(line[0:5])
        resname = line[5:10].strip() or "POL"
        atomname = line[10:15].strip() or "X"
        atomnr = int(line[15:20])
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        tail = line[44:].rstrip("\n")
    except Exception as exc:
        raise ValueError(f"Could not parse GRO atom line: {line!r}") from exc
    return GroAtom(resnr, resname, atomname, atomnr, np.array([x, y, z], dtype=float), tail)


def read_gro(gro_file: str | os.PathLike) -> tuple[str, list[GroAtom], np.ndarray]:
    """Read title, atoms and box lengths from a .gro file. Units are nm."""
    path = Path(gro_file)
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{gro_file} does not look like a valid GRO file.")

    title = lines[0]
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    if len(atom_lines) != n_atoms:
        raise ValueError(f"Expected {n_atoms} atoms in {gro_file}, found {len(atom_lines)}.")

    atoms = [_parse_gro_atom_line(line) for line in atom_lines]

    box_values = [float(x) for x in lines[2 + n_atoms].split()]
    if len(box_values) < 3:
        raise ValueError("GRO box line must contain at least three box lengths.")
    box_nm = np.array(box_values[:3], dtype=float)
    return title, atoms, box_nm


def write_gro(
    out_file: str | os.PathLike,
    atoms: Sequence[GroAtom],
    box_nm: Sequence[float],
    title: str = "BIO-SUSHY mini-bulk",
) -> str:
    """Write a simple orthorhombic .gro file. Units are nm."""
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    box = np.asarray(box_nm, dtype=float)
    with out_file.open("w") as handle:
        handle.write(f"{title}\n")
        handle.write(f"{len(atoms):5d}\n")
        for i, atom in enumerate(atoms, start=1):
            resnr = ((atom.resnr - 1) % 99999) + 1
            atomnr = ((i - 1) % 99999) + 1
            x, y, z = atom.xyz_nm
            handle.write(
                f"{resnr:5d}{atom.resname[:5]:>5s}{atom.atomname[:5]:>5s}"
                f"{atomnr:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
            )
        handle.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")
    return str(out_file)


def _random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """Generate a random 3D rotation matrix using a unit quaternion."""
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ],
        dtype=float,
    )


def _extract_single_molecule_mass_amu(top_file: str | os.PathLike) -> float | None:
    """
    Extract the total mass of one molecule from a Gromacs [ atoms ] section.

    Expected columns are usually:
    nr type resnr residue atom cgnr charge mass
    """
    lines = Path(top_file).read_text().splitlines()
    in_atoms = False
    total = 0.0
    found = False

    for raw in lines:
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[]").strip().lower()
            in_atoms = section == "atoms"
            continue
        if in_atoms:
            parts = line.split()
            if len(parts) >= 8:
                try:
                    total += float(parts[7])
                    found = True
                except ValueError:
                    pass

    return total if found and total > 0 else None


def _fallback_mass_from_gro_atoms(atoms: Sequence[GroAtom]) -> float:
    """Crude fallback mass estimate from atom names when topology masses are unavailable."""
    masses = {
        "H": 1.008,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.998,
        "P": 30.974,
        "S": 32.06,
        "CL": 35.45,
        "BR": 79.904,
        "I": 126.90,
    }
    total = 0.0
    for atom in atoms:
        name = "".join(ch for ch in atom.atomname.upper() if ch.isalpha())
        if name.startswith("CL"):
            key = "CL"
        elif name.startswith("BR"):
            key = "BR"
        else:
            key = name[:1] or "C"
        total += masses.get(key, 12.011)
    return total


def _update_molecules_count(top_file: str | os.PathLike, out_top: str | os.PathLike, n_chains: int) -> str:
    """
    Copy a topology file and update the first molecule count in [ molecules ].

    This assumes the single-chain topology contains one polymer molecule type.
    """
    lines = Path(top_file).read_text().splitlines()
    out_lines: list[str] = []
    in_molecules = False
    updated = False

    for raw in lines:
        stripped = raw.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped.strip("[]").strip().lower()
            in_molecules = section == "molecules"
            out_lines.append(raw)
            continue

        if in_molecules and not updated:
            data_part = raw.split(";", 1)[0].strip()
            if data_part:
                parts = data_part.split()
                if len(parts) >= 2:
                    comment = ""
                    if ";" in raw:
                        comment = " ;" + raw.split(";", 1)[1]
                    out_lines.append(f"{parts[0]:<18s}{int(n_chains):>6d}{comment}")
                    updated = True
                    continue
        out_lines.append(raw)

    if not updated:
        raise ValueError(
            "Could not find a molecule entry in the [ molecules ] section of the topology."
        )

    out_top = Path(out_top)
    out_top.parent.mkdir(parents=True, exist_ok=True)
    out_top.write_text("\n".join(out_lines) + "\n")
    return str(out_top)


def _copy_local_includes(top_file: str | os.PathLike, output_dir: str | os.PathLike) -> None:
    """Copy local include files referenced by a topology, when they exist."""
    top_path = Path(top_file)
    output_dir = Path(output_dir)
    for raw in top_path.read_text().splitlines():
        stripped = raw.strip()
        if stripped.startswith("#include"):
            token = stripped.replace("#include", "", 1).strip().strip('"<>')
            src = top_path.parent / token
            if src.exists() and src.is_file():
                dst = output_dir / token
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)


def build_polymer_bulk(
    single_chain_gro: str,
    single_chain_top: str,
    n_chains: int = 4,
    initial_density_g_cm3: float = 0.15,
    output_dir: str = "Bulk",
    random_seed: int = 42,
    padding_nm: float = 0.5,
) -> tuple[str, str]:
    """
    Build a small demonstrative periodic mini-bulk from one chain.

    Parameters
    ----------
    single_chain_gro, single_chain_top
        Gromacs coordinate and topology files for one polymer chain.
    n_chains
        Number of copies to place in the periodic box.
    initial_density_g_cm3
        Target initial density. Keep this low for tutorial robustness; the NPT
        protocol will compress the box later.
    output_dir
        Directory where bulk.gro and topol.top will be written.
    random_seed
        Reproducibility seed for rotations and small placement jitters.
    padding_nm
        Extra padding used when the density-derived box would be too small for
        the chain dimensions.

    Returns
    -------
    bulk_gro, bulk_top
    """
    if n_chains < 1:
        raise ValueError("n_chains must be at least 1.")
    if initial_density_g_cm3 <= 0:
        raise ValueError("initial_density_g_cm3 must be positive.")

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    title, base_atoms, _ = read_gro(single_chain_gro)
    base_xyz = np.array([atom.xyz_nm for atom in base_atoms], dtype=float)
    center = base_xyz.mean(axis=0)
    centered_xyz = base_xyz - center
    max_radius = float(np.max(np.linalg.norm(centered_xyz, axis=1)))

    mass_amu = _extract_single_molecule_mass_amu(single_chain_top)
    if mass_amu is None:
        mass_amu = _fallback_mass_from_gro_atoms(base_atoms)

    target_volume_nm3 = (mass_amu * n_chains * AMU_TO_G_PER_NM3_AT_DENSITY_1) / initial_density_g_cm3
    density_box_length = target_volume_nm3 ** (1.0 / 3.0)

    grid = math.ceil(n_chains ** (1.0 / 3.0))
    min_spacing = max(2.0 * max_radius + padding_nm, 0.8)
    box_length_nm = max(density_box_length, grid * min_spacing)
    spacing = box_length_nm / grid

    rng = np.random.default_rng(random_seed)
    grid_centers = []
    offset = spacing / 2.0
    for ix in range(grid):
        for iy in range(grid):
            for iz in range(grid):
                grid_centers.append(np.array([ix, iy, iz], dtype=float) * spacing + offset)
    rng.shuffle(grid_centers)

    all_atoms: list[GroAtom] = []
    atom_counter = 1
    for chain_idx in range(n_chains):
        rot = _random_rotation_matrix(rng)
        jitter = rng.uniform(-0.10 * spacing, 0.10 * spacing, size=3)
        chain_center = grid_centers[chain_idx] + jitter
        xyz = centered_xyz @ rot.T + chain_center
        xyz = np.mod(xyz, box_length_nm)

        for atom, coord in zip(base_atoms, xyz):
            all_atoms.append(
                GroAtom(
                    resnr=chain_idx + 1,
                    resname=atom.resname,
                    atomname=atom.atomname,
                    atomnr=atom_counter,
                    xyz_nm=coord,
                    tail=atom.tail,
                )
            )
            atom_counter += 1

    bulk_gro = write_gro(
        output / "bulk.gro",
        all_atoms,
        [box_length_nm, box_length_nm, box_length_nm],
        title=f"BIO-SUSHY mini-bulk: {n_chains} chains from {Path(single_chain_gro).name}",
    )
    bulk_top = _update_molecules_count(single_chain_top, output / "topol.top", n_chains)
    _copy_local_includes(single_chain_top, output)

    actual_initial_density = (mass_amu * n_chains * AMU_TO_G_PER_NM3_AT_DENSITY_1) / (box_length_nm**3)

    print("✅ Mini-bulk structure created.")
    print(f"   Chains:              {n_chains}")
    print(f"   Atoms:               {len(all_atoms)}")
    print(f"   Box length:          {box_length_nm:.3f} nm")
    print(f"   Initial density:     {actual_initial_density:.3f} g/cm^3")
    print(f"   Coordinates:         {bulk_gro}")
    print(f"   Topology:            {bulk_top}")

    return bulk_gro, bulk_top


def inspect_bulk_gro(gro_file: str, n_lines: int = 8) -> None:
    """Print a compact preview of a generated bulk .gro file."""
    title, atoms, box = read_gro(gro_file)
    print("\n" + "=" * 60)
    print("👀 INSPECTING GENERATED MINI-BULK")
    print("=" * 60)
    print(f"Title: {title}")
    print(f"Atoms: {len(atoms)}")
    print(f"Box:   {box[0]:.3f} × {box[1]:.3f} × {box[2]:.3f} nm")
    print("Preview:")
    for atom in atoms[:n_lines]:
        x, y, z = atom.xyz_nm
        print(f"   res {atom.resnr:>3d} {atom.resname:<5s} {atom.atomname:<5s} "
              f"{x:8.3f} {y:8.3f} {z:8.3f}")
    if len(atoms) > n_lines:
        print("   ...")
    print("-" * 60)
