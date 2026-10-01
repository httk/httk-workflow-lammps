"""Building blocks for the collect hooks of LAMMPS workflows."""

import re
from pathlib import Path

from httk.core import DataRecord
from httk.core.datastream.compression import open_compressed, split_compression_suffix

from .outputs import _ENERGY_FACTORS, _echo_state, average_total_energy_ev, parse_lammps_log

__all__ = ["find_outputs", "read_average_total_energy", "read_setup", "units_convertible"]

_AVERAGE_DEFINITION = "https://schemas.httk.org/defs/v0.1/properties/core/average_total_energy"
_AVERAGE_NAME = "_httk_average_total_energy"
_BANNER = re.compile(r"^LAMMPS \(")
_BANNER_LINES = 5
_ECHO_LINES = 2000


def _head(path: Path, limit: int) -> list[str]:
    """The first *limit* lines of *path*, decompressed; empty when it cannot be read."""

    lines: list[str] = []
    try:
        with path.open("rb") as raw, open_compressed(raw, compression="extension", name=path.name) as stream:
            for number, line in enumerate(stream):
                if number >= limit:
                    break
                lines.append(line.decode("utf-8", "replace"))
    except (OSError, EOFError):
        pass
    return lines


def find_outputs(directory: Path) -> tuple[Path, ...]:
    """Find the LAMMPS logs of a directory, sorted by name.

    Any file counts (``log.lammps`` has no ``.log`` suffix, and a log may be compressed)
    when a ``LAMMPS (`` banner line appears within its first five lines; only that head is read.
    LAMMPS prints the banner on the screen too, so when several files qualify, those that
    echo the input (the log does, a captured screen copy does not) are kept, then, if
    still several, those named ``log.*`` or ``*.log``.

    :param directory: The directory to look in.
    :return: The logs, sorted by name.
    """

    found = [
        path
        for path in sorted(Path(directory).iterdir())
        if path.is_file() and any(_BANNER.match(line) for line in _head(path, _BANNER_LINES))
    ]
    if len(found) > 1:
        found = [path for path in found if read_setup(path)[0] is not None] or found
    if len(found) > 1:
        named = [path for path in found if path.name.lower().startswith("log.") or _stem(path).endswith(".log")]
        found = named or found
    return tuple(found)


def _stem(path: Path) -> str:
    return split_compression_suffix(path.name)[0].lower()


def read_setup(path: Path) -> tuple[str | None, bool]:
    """Read the ``units`` style and the ``thermo_modify norm`` state echoed near the top of a LAMMPS log.

    The state is the last one in force within the first 2000 lines (a ``clear`` resets it); a later
    change is caught when the full parse of the collect step picks the state of the table itself.

    :param path: The LAMMPS log, optionally compressed.
    :return: The style (``None`` when the log echoes none, as with ``echo none``) and whether
        ``thermo_modify norm yes`` is echoed.
    """

    units, normalized = None, False
    for line in _head(path, _ECHO_LINES):
        units, normalized = _echo_state(units, normalized, line.split())
    return units, normalized


def units_convertible(units: str | None) -> bool:
    """Whether energies in the LAMMPS *units* style convert to eV.

    :param units: The style from :func:`read_setup`.
    :return: ``True`` for ``metal`` and ``real``.
    """

    return units in _ENERGY_FACTORS


def read_average_total_energy(path: Path) -> DataRecord:
    """Read the mean printed total energy of the last ``run``, in eV, from one LAMMPS log.

    LAMMPS prints no average, so this is the arithmetic mean of the ``TotEng`` column of the
    last ``run`` thermo table (``minimize`` tables excluded).

    :param path: The LAMMPS log, optionally compressed.
    :return: The mean as an ``average_total_energy`` property record.
    :raises ValueError: If the run did not complete, has no ``run`` thermo table, prints no
        ``TotEng`` column, uses ``units`` other than ``metal`` or ``real``, or has ``thermo_modify norm yes``.
    """

    result = parse_lammps_log(path)
    if not result.completed:
        raise ValueError(f"{path} is not a completed LAMMPS run")
    average = average_total_energy_ev(result)
    if average is None:
        raise ValueError(f"{path} has no run thermo table with TotEng; print etotal in thermo_style")
    return DataRecord.from_value(_AVERAGE_DEFINITION, _AVERAGE_NAME, average)
