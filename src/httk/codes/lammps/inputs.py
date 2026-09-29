"""Write LAMMPS data files from *httk* structures."""

import math
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

__all__ = ["write_lammps_data"]


def write_lammps_data(
    path: str | os.PathLike[str],
    *,
    structure: object,
    atom_style: str = "atomic",
    masses: Mapping[str, float] | None = None,
) -> Path:
    """Write a LAMMPS data file (for ``read_data``) of one periodic structure.

    The cell is rotated into the LAMMPS restricted-triclinic frame: ``a`` along
    *x*, ``b`` in the *xy* plane, and the box written as ``xlo xhi``, ``ylo
    yhi``, ``zlo zhi`` from the origin plus an ``xy xz yz`` tilt line unless the
    cell is orthogonal. Positions are Cartesian in Å in that frame. Atom types
    are numbered from 1 in order of first appearance among the sites, and each
    ``Masses`` line names its species in a comment. A data file is inherently
    floating point: a POSCAR/CONTCAR path is read as the floats it writes (its
    species are its element symbols), any other structure as floats of its exact
    values, and lengths are written rounded to 1e-12 Å. This needs
    *httk-atomistic* (the ``atomistic`` extra).

    :param path: Write the data file to this path.
    :param structure: The structure: a POSCAR/CIF path or anything
        ``httk.atomistic.UnitcellStructureView`` accepts.
    :param atom_style: The LAMMPS atom style of the ``Atoms`` section; only ``atomic`` is supported.
    :param masses: Map species names to masses (in the script's mass unit, e.g.
        g/mol for ``units metal``). Without it the structure's own species masses
        are used; when some species has neither, no ``Masses`` section is written
        and the script must set them with ``mass``.
    :return: The written path.
    :raises ValueError: If the atom style is unsupported, a given mass table misses a
        species, the cell is left-handed, or a POSCAR has no species line.
    """

    # ponytail: only `atomic`; add charge/molecular columns when a workflow needs them.
    if atom_style != "atomic":
        raise ValueError(f"atom_style must be 'atomic', not {atom_style!r}")
    cell, sites, positions, own = _read_structure(structure)
    a, b, c = cell
    if _dot(a, _cross(b, c)) <= 0:
        raise ValueError("the cell is left-handed; LAMMPS needs a right-handed cell (reorder or negate a vector)")
    # The standard restricted-triclinic frame (LAMMPS Howto triclinic): x along a, y in the ab plane.
    x = _unit(a)
    y = _unit([bi - _dot(b, x) * xi for bi, xi in zip(b, x, strict=True)])
    z = _cross(x, y)
    lx, xy, ly, xz, yz, lz = _dot(a, x), _dot(b, x), _dot(b, y), _dot(c, x), _dot(c, y), _dot(c, z)

    names = tuple(dict.fromkeys(sites))
    types = {name: index for index, name in enumerate(names, start=1)}
    if masses is not None:
        missing = [name for name in names if name not in masses]
        if missing:
            raise ValueError(f"no mass given for species {', '.join(missing)}")
        table: Mapping[str, float] | None = masses
    else:
        table = own if all(name in own for name in names) else None

    lines = [
        f"LAMMPS data file written by httk-workflow-lammps; types: {' '.join(names)}",
        "",
        f"{len(sites)} atoms",
        f"{len(names)} atom types",
        "",
        f"0.0 {_f(lx)} xlo xhi",
        f"0.0 {_f(ly)} ylo yhi",
        f"0.0 {_f(lz)} zlo zhi",
    ]
    if _f(xy) != "0.0" or _f(xz) != "0.0" or _f(yz) != "0.0":
        lines.append(f"{_f(xy)} {_f(xz)} {_f(yz)} xy xz yz")
    if table is not None:
        lines += ["", "Masses", ""] + [f"{types[name]} {float(table[name])!r} # {name}" for name in names]
    lines += ["", f"Atoms # {atom_style}", ""]
    for index, (name, r) in enumerate(zip(sites, positions, strict=True), start=1):
        lines.append(f"{index} {types[name]} {_f(_dot(r, x))} {_f(_dot(r, y))} {_f(_dot(r, z))}")
    destination = Path(path)
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination


def _read_structure(
    structure: object,
) -> tuple[list[list[float]], tuple[str, ...], list[list[float]], dict[str, float]]:
    """Return the cell rows, site species names, Cartesian positions (Å) and known species masses.

    A POSCAR/CONTCAR path is read as the floats it writes: loading it as an exact
    structure would turn each short decimal into the simplest fraction within its
    written precision (``-1.92`` into ``-23/12``).
    """

    import httk.core
    from httk.atomistic import UnitcellStructureView  # pyright: ignore[reportMissingImports]
    from httk.core.loading import adapt_result

    if isinstance(structure, str | os.PathLike):
        # The precision only answers the POSCAR reader's symmetry-tolerance advice; the
        # structure readers take it uniformly and nothing below uses it.
        data = httk.core.load(os.fspath(structure), raw=True, precision=1e-6)
        if isinstance(data, Mapping) and data.get("format") == "vasp-poscar":
            return _poscar_floats(data)
        structure = adapt_result(data, raw=False)
    view = UnitcellStructureView(cast(Any, structure))
    cell = [[float(x) for x in row] for row in view.lattice_vectors]
    positions = [[float(x) for x in row] for row in view.cartesian_site_positions]
    own = {item.name: float(item.mass[0]) for item in view.species if item.mass}
    return cell, tuple(view.species_at_sites), positions, own


def _poscar_floats(
    data: Mapping[str, Any],
) -> tuple[list[list[float]], tuple[str, ...], list[list[float]], dict[str, float]]:
    if data["symbols"] is None:
        raise ValueError("a VASP-4 POSCAR has no species line; give a VASP-5 file")
    raw = [[float(x) for x in row] for row in data["cell"]]
    if data["scale"] is not None:
        scale = float(data["scale"])
    elif data["volume"] is not None:
        scale = (float(data["volume"]) / abs(_dot(raw[0], _cross(raw[1], raw[2])))) ** (1 / 3)
    else:
        scale = 1.0
    cell = [[scale * x for x in row] for row in raw]
    coords = [[float(x) for x in row] for row in data["coords"]]
    if data["cartesian"]:
        positions = [[scale * x for x in row] for row in coords]
    else:
        positions = [[sum(f * cell[i][k] for i, f in enumerate(row)) for k in range(3)] for row in coords]
    sites = tuple(symbol for symbol, count in zip(data["symbols"], data["counts"], strict=True) for _ in range(count))
    return cell, sites, positions, {}


def _f(value: float) -> str:
    # Round away float noise below 1e-12 Å (and the sign of a negative zero).
    return repr(round(value, 12) + 0.0)


def _unit(u: list[float]) -> list[float]:
    norm = math.sqrt(_dot(u, u))
    return [x / norm for x in u]


def _dot(u: list[float], v: list[float]) -> float:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def _cross(u: list[float], v: list[float]) -> list[float]:
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
