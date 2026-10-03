"""Recognize hook for ``lammps.calculation``: one LAMMPS log with its input script beside it."""

from pathlib import Path

from httk.core.datastream.compression import split_compression_suffix
from httk.workflow.calculations import content_digest
from httk.workflow.collecting import existing_file
from httk.workflow.hookapi import Claim, Unclaimed

from httk.codes.lammps.collect import _head, find_outputs, read_setup, units_convertible
from httk.codes.lammps.outputs import NORMALIZED_REASON

_FIRST_ARGUMENT = ("read_data", "read_restart", "include")
_ALL_ARGUMENTS = ("pair_coeff", "molecule")


def _referenced(directory: Path, script: str) -> list[str]:
    """The plain relative files in *directory* that the script names.

    These are the first argument of ``read_data``/``read_restart``/``include`` and any argument
    of ``pair_coeff`` (potential files) and ``molecule`` lines.
    """

    names: list[str] = []
    path = existing_file(directory / script)
    for line in _head(path, 1_000_000) if path else ():
        fields = line.split("#")[0].split()
        if fields[:1] and fields[0] in (*_FIRST_ARGUMENT, *_ALL_ARGUMENTS):
            for name in fields[1:2] if fields[0] in _FIRST_ARGUMENT else fields[1:]:
                name = name.strip("\"'")
                if "$" not in name and "/" not in name and existing_file(directory / name):
                    names.append(name)
    return names


def recognize(directory: Path) -> Claim | Unclaimed | None:
    """Claim a directory holding exactly one LAMMPS log with convertible units and one input script.

    :param directory: The directory to examine.
    :return: A claim identified by the input script's content and the files it reads, ``Unclaimed`` for a LAMMPS
        directory that cannot be collected, or ``None`` when it holds no LAMMPS log.
    """

    logs = list(find_outputs(directory))
    if not logs:
        return None
    if len(logs) > 1:
        return Unclaimed(f"several LAMMPS logs: {', '.join(path.name for path in logs)}")
    units, normalized = read_setup(logs[0])
    if not units_convertible(units):
        return Unclaimed(f"LAMMPS units {units or 'unknown'} cannot be converted to eV")
    if normalized is None:
        return Unclaimed("LAMMPS thermo normalization state is unknown")
    if normalized:
        return Unclaimed(NORMALIZED_REASON)
    stripped = (split_compression_suffix(path.name)[0] for path in directory.iterdir() if path.is_file())
    scripts = sorted(
        name
        for name in stripped
        if name != split_compression_suffix(logs[0].name)[0]
        and (name.lower().startswith("in.") or name.lower().endswith(".in"))
    )
    if len(scripts) != 1:
        return Unclaimed(
            f"no in.* or *.in input script beside {logs[0].name}"
            if not scripts
            else f"several LAMMPS input scripts: {', '.join(scripts)}"
        )
    return Claim(content_digest(directory, [*scripts, *_referenced(directory, scripts[0])]))
