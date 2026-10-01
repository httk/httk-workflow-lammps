"""Parse the log of a LAMMPS (``lmp``) run.

Pure stdlib parsing: nothing here runs a program or imports *httk* code, so a
result can be read anywhere the log file is.
"""

import os
import re
import statistics
from dataclasses import dataclass, replace
from pathlib import Path

from httk.core.datastream.compression import open_compressed

__all__ = ["KCAL_MOL_TO_EV", "LammpsResult", "ThermoTable", "average_total_energy_ev", "parse_lammps_log"]

#: One kcal/mol in electronvolts: the thermochemical kcal (4.184 kJ) times the
#: CODATA kJ/mol to eV factor, the same basis as ``KJ_MOL_TO_EV`` of httk-workflow-gromacs.
KCAL_MOL_TO_EV: float = 4.184 * 0.010364269656262175
NORMALIZED_REASON = "thermo_modify norm yes reports per-atom energies; the average total energy is for the whole system"
_ENERGY_FACTORS = {"metal": 1.0, "real": KCAL_MOL_TO_EV}
_LOOP_STEPS = re.compile(r"^Loop time of .* for (\d+) steps")

_STOP = "Stopping criterion ="
#: The minimization stopping criteria that mean the tolerances were not reached.
_UNCONVERGED_STOPS = frozenset({"max iterations", "max force evaluations"})


@dataclass(frozen=True)
class ThermoTable:
    """The thermodynamic output of one ``run`` or ``minimize`` command.

    :param columns: The column names of the header line, e.g. ``("Step", "Temp", "PotEng")``.
    :param rows: One tuple of values per printed step, as floats.
    :param kind: ``"run"`` or ``"minimize"``, from the echoed command before the table or the
        ``Minimization stats:`` block after it.
    :param steps: The step count of the ``Loop time of`` line, or ``None`` for a table cut short.
    :param units: The ``units`` style in force when the table was printed (the last echoed ``units``
        command since any ``clear``), or ``None`` when none was echoed.
    :param normalized: Whether ``thermo_modify norm yes`` was in force (energies per atom).
    """

    columns: tuple[str, ...]
    rows: tuple[tuple[float, ...], ...]
    kind: str = "run"
    steps: int | None = None
    units: str | None = None
    normalized: bool = False

    @property
    def last(self) -> dict[str, float]:
        """The last row as a column-name mapping, empty when the table has no rows."""
        return dict(zip(self.columns, self.rows[-1], strict=True)) if self.rows else {}


@dataclass(frozen=True)
class LammpsResult:
    """What one LAMMPS log says about its run.

    :param thermo: The thermodynamic output tables, one per ``run``/``minimize`` command, in order.
    :param completed: Whether LAMMPS reached its normal end (the ``Total wall time:`` line).
    :param errors: The ``ERROR:`` and ``ERROR on proc N:`` lines, in order and without repeats.
    :param warnings: The ``WARNING:`` lines, in order and without repeats.
    :param minimization_stop: The stopping criterion of the last minimization, or ``None``.
    :param units: The last ``units`` style echoed in the log (``None`` after a ``clear``, or when none is echoed).
    :param normalized: Whether ``thermo_modify norm yes`` is in force at the end of the log.
    """

    thermo: tuple[ThermoTable, ...]
    completed: bool
    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    minimization_stop: str | None
    units: str | None = None
    normalized: bool = False

    @property
    def minimization_converged(self) -> bool | None:
        """Whether the last minimization met its tolerances, or ``None`` when there was none.

        LAMMPS names no convergence verdict; stopping on ``max iterations`` or
        ``max force evaluations`` counts as not converged, any other criterion as converged.
        """
        return None if self.minimization_stop is None else self.minimization_stop not in _UNCONVERGED_STOPS


def parse_lammps_log(path: str | os.PathLike[str]) -> LammpsResult:
    """Parse one LAMMPS log file (``-log``, default ``log.lammps``).

    :param path: Read the LAMMPS log saved at this path, optionally compressed.
    :return: The parsed result.
    :raises FileNotFoundError: If the log file does not exist.
    """

    path = Path(path)
    with path.open("rb") as raw, open_compressed(raw, compression="extension", name=path.name) as stream:
        return _parse(stream.read().decode("utf-8", errors="replace"))


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""


def _floats(fields: list[str]) -> tuple[float, ...] | None:
    try:
        return tuple(float(field) for field in fields)
    except ValueError:
        return None


def _echo_state(units: str | None, normalized: bool, fields: list[str]) -> tuple[str | None, bool]:
    """Update the (units, norm) state in force by one echoed input line, given as split fields."""

    match fields:
        case ["units", style, *_]:
            return style, normalized
        case ["clear"]:
            return None, False
        case ["thermo_modify", *args] if "norm" in args[:-1]:
            last = max(i for i, arg in enumerate(args[:-1]) if arg == "norm")
            return units, args[last + 1].lower() == "yes"
    return units, normalized


def _parse(text: str, screen: str = "") -> LammpsResult:
    """Parse a log; *screen* (captured standard output) contributes only its error lines.

    An error raised before LAMMPS opens its log, such as an unreadable input
    script, is printed to the screen alone.
    """

    tables: list[ThermoTable] = []
    columns: tuple[str, ...] | None = None
    rows: list[tuple[float, ...]] = []
    stop: str | None = None
    command = "run"  # the kind of the last echoed run/minimize command
    units: str | None = None
    normalized = False
    # ponytail: `thermo_style multi`/`yaml` blocks are not tables here; parse them when a workflow needs them.
    for line in text.splitlines():
        fields = line.split()
        if columns is not None:
            values = _floats(fields) if len(fields) == len(columns) else None
            if values is not None:
                rows.append(values)
            elif line.startswith("Loop time of"):
                loop = _LOOP_STEPS.match(line)
                tables.append(
                    ThermoTable(columns, tuple(rows), command, int(loop.group(1)) if loop else None, units, normalized)
                )
                columns = None
            # anything else (a WARNING, a fix's message) interleaves with the rows and is skipped
        elif fields[:1] in (["units"], ["clear"], ["thermo_modify"]):
            units, normalized = _echo_state(units, normalized, fields)
        elif fields[:1] == ["Step"]:
            columns, rows = tuple(fields), []
        elif line.strip().startswith(_STOP):
            stop = line.split("=", 1)[1].strip()
        elif fields[:1] in (["run"], ["minimize"]):
            command = fields[0]
        elif line.startswith("Minimization stats:") and tables:
            tables[-1] = replace(tables[-1], kind="minimize")
    if columns is not None:  # a run cut short keeps the rows it printed
        tables.append(ThermoTable(columns, tuple(rows), command, None, units, normalized))
    lines = [line.strip() for line in (text + "\n" + screen).splitlines()]
    return LammpsResult(
        thermo=tuple(tables),
        completed="Total wall time:" in text,
        errors=tuple(dict.fromkeys(line for line in lines if line.startswith(("ERROR:", "ERROR on proc")))),
        warnings=tuple(dict.fromkeys(line for line in lines if line.startswith("WARNING:"))),
        minimization_stop=stop,
        units=units,
        normalized=normalized,
    )


def average_total_energy_ev(result: LammpsResult) -> float | None:
    """Mean of the printed ``TotEng`` over the last ``run`` thermo table, in eV.

    LAMMPS prints no average itself; this is the arithmetic mean of the printed rows.
    ``minimize`` tables and zero-step ``run 0`` tables are skipped.

    :param result: The parsed log.
    :return: The mean in eV, or ``None`` when there is no ``run`` table or it has no ``TotEng`` rows.
    :raises ValueError: If the echoed ``units`` style is not ``metal`` or ``real``, or is not echoed,
        or the table's input state has ``thermo_modify norm yes``. The units and norm state are those in force
        for the chosen table.
    """

    runs = [table for table in result.thermo if table.kind == "run" and table.steps != 0]
    if not runs or "TotEng" not in runs[-1].columns or not runs[-1].rows:
        return None
    table = runs[-1]
    if table.normalized:
        raise ValueError(NORMALIZED_REASON)
    if table.units not in _ENERGY_FACTORS:
        raise ValueError(f"LAMMPS units {table.units or 'unknown (no echoed units line)'!r} cannot be converted to eV")
    index = table.columns.index("TotEng")
    return statistics.fmean(row[index] for row in table.rows) * _ENERGY_FACTORS[table.units]
