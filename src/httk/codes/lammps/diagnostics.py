"""Classify a finished LAMMPS run into stable diagnostics."""

import os
from pathlib import Path

from httk.workflow.codes import Diagnostic

from .outputs import _parse, _read

__all__ = ["diagnose_lammps"]


def diagnose_lammps(
    directory: str | os.PathLike[str] = ".", *, log: str = "log.lammps", output: str = "lammps.out"
) -> tuple[Diagnostic, ...]:
    """Diagnose a LAMMPS run from its log and captured standard output.

    The codes are stable: ``lammps.error`` (fatal; LAMMPS stopped with an
    ``ERROR`` line, which is the summary), ``lammps.incomplete`` (error; no
    ``Total wall time:`` line and no error, e.g. a killed process) and
    ``lammps.minimization_not_converged`` (error; the last minimization stopped
    on ``max iterations`` or ``max force evaluations``). ``WARNING`` lines are
    not diagnostics; they are kept in :attr:`LammpsResult.warnings
    <httk.codes.lammps.LammpsResult.warnings>`. A clean run has none. A missing
    file is diagnosed like an empty one.

    :param directory: Read the run files from this directory.
    :param log: The name of the LAMMPS log in *directory*.
    :param output: The name of the saved standard output in *directory*; only its errors are read.
    :return: The diagnostics, empty for a clean run.
    """

    root = Path(directory)
    result = _parse(_read(root / log), _read(root / output))
    diagnostics: list[Diagnostic] = []
    if result.errors:
        diagnostics.append(Diagnostic("lammps.error", "fatal", result.errors[0], log, "\n".join(result.errors)))
    elif not result.completed:
        diagnostics.append(Diagnostic("lammps.incomplete", "error", f"{log} has no Total wall time: line", log))
    if result.minimization_converged is False:
        diagnostics.append(
            Diagnostic(
                "lammps.minimization_not_converged",
                "error",
                f"minimization stopped on {result.minimization_stop}",
                log,
            )
        )
    return tuple(diagnostics)
