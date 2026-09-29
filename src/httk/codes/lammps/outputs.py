"""Parse the log of a LAMMPS (``lmp``) run.

Pure stdlib parsing: nothing here runs a program or imports *httk* code, so a
result can be read anywhere the log file is.
"""

import os
from dataclasses import dataclass
from pathlib import Path

__all__ = ["LammpsResult", "ThermoTable", "parse_lammps_log"]

_STOP = "Stopping criterion ="
#: The minimization stopping criteria that mean the tolerances were not reached.
_UNCONVERGED_STOPS = frozenset({"max iterations", "max force evaluations"})


@dataclass(frozen=True)
class ThermoTable:
    """The thermodynamic output of one ``run`` or ``minimize`` command.

    :param columns: The column names of the header line, e.g. ``("Step", "Temp", "PotEng")``.
    :param rows: One tuple of values per printed step, as floats.
    """

    columns: tuple[str, ...]
    rows: tuple[tuple[float, ...], ...]

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
    """

    thermo: tuple[ThermoTable, ...]
    completed: bool
    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    minimization_stop: str | None

    @property
    def minimization_converged(self) -> bool | None:
        """Whether the last minimization met its tolerances, or ``None`` when there was none.

        LAMMPS names no convergence verdict; stopping on ``max iterations`` or
        ``max force evaluations`` counts as not converged, any other criterion as converged.
        """
        return None if self.minimization_stop is None else self.minimization_stop not in _UNCONVERGED_STOPS


def parse_lammps_log(path: str | os.PathLike[str]) -> LammpsResult:
    """Parse one LAMMPS log file (``-log``, default ``log.lammps``).

    :param path: Read the LAMMPS log saved at this path.
    :return: The parsed result.
    :raises FileNotFoundError: If the log file does not exist.
    """

    return _parse(Path(path).read_text(encoding="utf-8", errors="replace"))


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""


def _floats(fields: list[str]) -> tuple[float, ...] | None:
    try:
        return tuple(float(field) for field in fields)
    except ValueError:
        return None


def _parse(text: str, screen: str = "") -> LammpsResult:
    """Parse a log; *screen* (captured standard output) contributes only its error lines.

    An error raised before LAMMPS opens its log, such as an unreadable input
    script, is printed to the screen alone.
    """

    tables: list[ThermoTable] = []
    columns: tuple[str, ...] | None = None
    rows: list[tuple[float, ...]] = []
    stop: str | None = None
    # ponytail: `thermo_style multi`/`yaml` blocks are not tables here; parse them when a workflow needs them.
    for line in text.splitlines():
        fields = line.split()
        if columns is not None:
            values = _floats(fields) if len(fields) == len(columns) else None
            if values is not None:
                rows.append(values)
            elif line.startswith("Loop time of"):
                tables.append(ThermoTable(columns, tuple(rows)))
                columns = None
            # anything else (a WARNING, a fix's message) interleaves with the rows and is skipped
        elif fields[:1] == ["Step"]:
            columns, rows = tuple(fields), []
        elif line.strip().startswith(_STOP):
            stop = line.split("=", 1)[1].strip()
    if columns is not None:  # a run cut short keeps the rows it printed
        tables.append(ThermoTable(columns, tuple(rows)))
    lines = [line.strip() for line in (text + "\n" + screen).splitlines()]
    return LammpsResult(
        thermo=tuple(tables),
        completed="Total wall time:" in text,
        errors=tuple(dict.fromkeys(line for line in lines if line.startswith(("ERROR:", "ERROR on proc")))),
        warnings=tuple(dict.fromkeys(line for line in lines if line.startswith("WARNING:"))),
        minimization_stop=stop,
    )
