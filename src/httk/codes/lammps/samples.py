"""Join explicitly selected LAMMPS dump segments and thermo tables.

The dump is read by the exact ``LammpsTrajectory`` backend of the optional *httk-atomistic*.
"""

import math
import os
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal, cast

from httk.core import definition_ids, load_property_definition
from httk.core.units import default_registry

from .outputs import parse_lammps_log

__all__ = ["LammpsSample", "ThermoConversion", "lammps_samples"]

_EV_J = 1.602176634e-19
_KB_EV = 8.617333262145e-5
_STRESS = ("Pxx", "Pyy", "Pzz", "Pyz", "Pxz", "Pxy")  # Voigt order
# No calorie/mole unit upstream: exact thermochemical kcal/mol (4184 J / N_A) in eV.
_KCAL_MOL_EV = 4184 / (6.02214076e23 * _EV_J)


@dataclass(frozen=True, slots=True)
class ThermoConversion:
    """Explicit meaning and conversion of a thermo column.

    :param definition: IRI of the OPTIMADE property definition the value is reported as.
    :param scale: Multiply the printed value by this finite scale to obtain the definition's unit; choosing it is the caller's responsibility.
    :param per_atom: If true, additionally multiply by the frame's atom count.
    """

    definition: str
    scale: float
    per_atom: bool = False

    def __post_init__(self) -> None:
        """Require a definition and a finite explicit scale."""
        if not self.definition or not math.isfinite(self.scale):
            raise ValueError("conversions require a definition and a finite scale")


@dataclass(frozen=True, slots=True)
class LammpsSample:
    """One joined structure and its converted thermodynamic values.

    :param structure: Exact ``httk.atomistic.UnitcellStructure`` from the selected dump segment.
    :param step: Integer step shared by the dump and selected table.
    :param thermo: Pairs of property definition name and value in the definition's unit, in caller column order.
    :param observables: Requested trajectory observables in caller order.
    """

    structure: Any
    step: int
    thermo: tuple[tuple[str, Any], ...]
    observables: tuple[Any, ...]


def lammps_samples(
    dump: str | os.PathLike[str],
    log: str | os.PathLike[str],
    *,
    species: Mapping[int, str],
    units: str,
    dimension: int,
    segment: int,
    table_index: int,
    columns: Sequence[str],
    observables: Sequence[str] = (),
    join: Literal["intersection", "strict"] = "intersection",
    total_atom_count: int | None = None,
    stress_basis_matches_dump: bool = False,
    custom_columns: Mapping[str, ThermoConversion] | None = None,
    trajectory_options: Mapping[str, Any] | None = None,
) -> Iterator[LammpsSample]:
    """Stream a step join without interpolation or implicit run selection.

    Requires *httk-atomistic* for the exact dump trajectory. Both segment and thermo table must
    be explicitly selected from the same physical run. An intersection emits
    shared steps only; a strict join requires identical ordered step sequences.
    A strict mismatch can be detected after earlier samples have been yielded.
    Log tables are materialized by this package's log parser; dump frames stream.

    :param dump: Native text dump path, optionally compressed.
    :param log: LAMMPS log path, optionally compressed.
    :param species: Explicit atom-type to element mapping.
    :param units: Explicit metal, real, or lj source units, checked against metadata.
    :param dimension: Explicit simulation dimension; this bulk converter supports only 3.
    :param segment: Nonnegative dump segment index.
    :param table_index: Nonnegative thermo table index.
    :param columns: Selected thermo column names in output order. The six stress components Pxx, Pyy, Pzz, Pyz, Pxz, Pxy must be selected together and yield one stress tensor value at the first one's position.
    :param observables: Additional canonical trajectory observable names to stream.
    :param join: Exact step intersection or identical-step strict join.
    :param total_atom_count: Global atom count for normalized energies when the log omits Atoms.
    :param stress_basis_matches_dump: Explicit assertion that requested stress components use the dump Cartesian basis.
    :param custom_columns: Explicit meanings/scales for nonstandard selected columns.
    :param trajectory_options: Additional LammpsTrajectory arguments, such as timestep or lj_scales.
    :yields: Joined immutable samples in dump order.
    :raises ImportError: If *httk-atomistic* is unavailable.
    :raises ValueError: If metadata, step identity, selection, units, or conversions are invalid, or stress components are incomplete.
    """
    try:
        from httk.atomistic.integrations.lammps.trajectory import (  # pyright: ignore[reportMissingImports]
            LammpsTrajectory,
        )
    except ImportError as exc:
        raise ImportError("LAMMPS dump/log joins require httk-atomistic") from exc
    if join not in ("intersection", "strict"):
        raise ValueError("join must be 'intersection' or 'strict'")
    if dimension != 3:
        raise ValueError("bulk thermo conversion requires dimension=3; 2D pressure and volume have different units")
    if any(isinstance(index, bool) or not isinstance(index, int) or index < 0 for index in (segment, table_index)):
        raise ValueError("segment and table_index must be nonnegative integers")
    if isinstance(columns, str) or isinstance(observables, str):
        raise ValueError("columns and observables must be sequences of names")
    if total_atom_count is not None and (
        isinstance(total_atom_count, bool) or not isinstance(total_atom_count, int) or total_atom_count <= 0
    ):
        raise ValueError("total_atom_count must be a positive integer")
    selected = tuple(columns)
    if not selected or len(set(selected)) != len(selected):
        raise ValueError("columns must be nonempty and distinct")
    result = parse_lammps_log(log)
    if table_index >= len(result.thermo):
        raise ValueError("table_index is outside the parsed thermo tables")
    table = result.thermo[table_index]
    if table.units != units:
        raise ValueError("thermo units must be known and agree with explicit dump units")
    if len(set(table.columns)) != len(table.columns) or "Step" not in table.columns:
        raise ValueError("thermo table requires distinct columns including Step")
    if any(name not in table.columns for name in selected):
        raise ValueError("a selected column is absent from the thermo table")
    options = dict(trajectory_options or {})
    if {"species", "units", "segment", "source"}.intersection(options):
        raise ValueError("trajectory_options cannot override explicit source, species, units, or segment")
    trajectory = LammpsTrajectory(dump, species=species, units=units, segment=segment, **options)
    factors = _factors(units, options.get("lj_scales"))
    custom = dict(custom_columns or {})
    if any(c.definition == definition_ids.STRESS_TENSOR for c in custom.values()):
        raise ValueError("the stress tensor comes only from the six Pxx..Pxy columns, not a custom conversion")
    stress = [name for name in selected if name in _STRESS and name not in custom]
    if stress:
        if missing := [name for name in _STRESS if name not in stress]:
            raise ValueError(f"stress tensor requires all six components; missing {missing}")
        if not stress_basis_matches_dump:
            raise ValueError("stress conversion requires explicit stress_basis_matches_dump=True")
    # One plan entry per output value: (conversion, column indices); a stress tensor has six.
    plan: list[tuple[ThermoConversion, tuple[int, ...]]] = []
    for name in selected:
        if name in stress:
            if name == stress[0]:
                tensor = ThermoConversion(definition_ids.STRESS_TENSOR, -factors[2])
                plan.append((tensor, tuple(table.columns.index(c) for c in _STRESS)))
        else:
            conversion = custom.get(name) or _standard(name, factors, table.normalized)
            plan.append((conversion, (table.columns.index(name),)))
    conversions = tuple(conversion for conversion, _ in plan)
    names = tuple(load_property_definition(conversion.definition).name for conversion in conversions)
    if len(set(names)) != len(names):
        raise ValueError("converted property names must be distinct")
    indices = tuple(table.columns.index(name) for name in selected)
    step_index = table.columns.index("Step")
    rows: dict[int, tuple[float, ...]] = {}
    for table_row in table.rows:
        value = table_row[step_index]
        if not math.isfinite(value) or not value.is_integer() or abs(value) >= 2**53:
            raise ValueError("thermo Step must be an exactly representable integer below 2**53")
        step = int(value)
        if rows and step <= next(reversed(rows)):
            raise ValueError("selected thermo table must have strictly increasing unique steps")
        if any(not math.isfinite(table_row[index]) for index in indices):
            raise ValueError("selected thermo values must be finite")
        rows[step] = table_row
    expected = iter(rows)
    for frame, values in trajectory.samples("step", *observables):
        step = values[0]
        if join == "strict" and next(expected, None) != step:
            raise ValueError("strict dump/thermo step sequences differ")
        row = rows.get(step)
        if row is None:
            continue
        atom_count = total_atom_count
        if "Atoms" in table.columns:
            printed_count = row[table.columns.index("Atoms")]
            if not math.isfinite(printed_count) or not printed_count.is_integer() or printed_count <= 0:
                raise ValueError("thermo Atoms must be a positive integer")
            if atom_count is not None and atom_count != int(printed_count):
                raise ValueError("total_atom_count disagrees with thermo Atoms")
            atom_count = int(printed_count)
        if any(conversion.per_atom for conversion in conversions) and (
            atom_count is None or atom_count < len(trajectory.species_at_sites)
        ):
            raise ValueError("normalized energy requires a global atom count, from Atoms or total_atom_count")
        converted = tuple(
            (name, _scaled(conversion, [row[index] for index in group], cast(int, atom_count)))
            for name, (conversion, group) in zip(names, plan, strict=True)
        )
        flat = [x for _, value in converted for x in (value if isinstance(value, tuple) else (value,))]
        if not all(math.isfinite(x) for x in flat):
            raise ValueError("converted thermo values are nonfinite")
        yield LammpsSample(frame, step, converted, values[1:])
    if join == "strict" and next(expected, None) is not None:
        raise ValueError("strict dump/thermo step sequences differ")


def _scaled(conversion: ThermoConversion, printed: list[float], atom_count: int) -> float | tuple[float, ...]:
    scale = conversion.scale * (atom_count if conversion.per_atom else 1)
    values = [value * scale for value in printed]
    return tuple(values) if conversion.definition == definition_ids.STRESS_TENSOR else values[0]


def _factors(units: str, lj_scales: Any) -> tuple[float, float, float, float]:
    """Return length, energy, pressure and temperature factors to angstrom, eV, GPa and K."""
    factor = default_registry().factor
    if units == "metal":
        return 1.0, 1.0, float(factor("bar", "GPa").factor), 1.0
    if units == "real":
        return 1.0, _KCAL_MOL_EV, float(factor("atm", "GPa").factor), 1.0
    length, _, energy = (float(value) for value in lj_scales)
    return length, energy, energy / length**3 * float(factor("angstrom^-3*eV", "GPa").factor), energy / _KB_EV


def _standard(name: str, factors: tuple[float, float, float, float], normalized: bool | None) -> ThermoConversion:
    length, energy, pressure, temperature = factors
    energies = {
        "TotEng": definition_ids.TOTAL_ENERGY,
        "PotEng": definition_ids.POTENTIAL_ENERGY,
        "KinEng": definition_ids.KINETIC_ENERGY,
        "Enthalpy": definition_ids.ENTHALPY,
    }
    if name in energies:
        if normalized is None:
            raise ValueError("energy normalization is unknown")
        return ThermoConversion(energies[name], energy, normalized)
    if name == "Temp":
        return ThermoConversion(definition_ids.TEMPERATURE, temperature)
    if name == "Volume":
        return ThermoConversion(definition_ids.VOLUME, length**3)
    if name == "Press":
        return ThermoConversion(definition_ids.PRESSURE, pressure)
    raise ValueError(f"column {name!r} requires an explicit ThermoConversion")
