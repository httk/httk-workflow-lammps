"""Supervised LAMMPS execution and its classified run report."""

import dataclasses
import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from httk.workflow.codes import Diagnostic, ProcessReport, ProcessSupervisor, launch_command, write_json_atomic

from .diagnostics import diagnose_lammps
from .outputs import LammpsResult, _parse, _read

__all__ = ["LammpsRunReport", "run_lammps"]


@dataclass(frozen=True)
class LammpsRunReport:
    """Classified result of one supervised LAMMPS execution.

    The classification is one of ``completed``, ``crashed`` (a ``lammps.error``
    diagnostic), ``nonconverged`` (a minimization that did not meet its
    tolerances), ``process_failure`` (a nonzero exit or an incomplete log) and
    ``timeout``.

    :param process: The supervised process result.
    :param classification: The final run classification.
    :param diagnostics: The diagnostics of the finished run.
    :param result: What the LAMMPS log says.
    """

    process: ProcessReport
    classification: str
    diagnostics: tuple[Diagnostic, ...]
    result: LammpsResult

    @property
    def ok(self) -> bool:
        """Whether the run completed cleanly."""
        return self.classification == "completed"

    def as_mapping(self) -> dict[str, object]:
        """Serialize the report for JSON storage.

        :return: The JSON-compatible report mapping.
        """
        return {
            "format": "httk-lammps-run-report",
            "format_version": 1,
            "process": self.process.as_mapping(),
            "classification": self.classification,
            "diagnostics": [item.as_mapping() for item in self.diagnostics],
            "result": dataclasses.asdict(self.result),
        }

    def write(self, path: str | os.PathLike[str]) -> Path:
        """Write the report as JSON.

        :param path: Write the report to this path.
        :return: The report path.
        """
        destination = Path(path)
        write_json_atomic(destination, self.as_mapping())
        return destination


def run_lammps(
    argv: Sequence[str],
    *,
    directory: str | os.PathLike[str] = ".",
    input_file: str = "in.lammps",
    log_file: str = "log.lammps",
    output_file: str = "lammps.out",
    timeout: float | None = None,
    launch: bool | None = None,
    termination_grace: float = 10.0,
    report_path: str | os.PathLike[str] = "lammps-run-report.json",
) -> LammpsRunReport:
    """Run LAMMPS under supervision and write a classified report.

    *argv* names the program (for example ``["lmp"]``); ``-in INPUT_FILE -log LOG_FILE``
    is appended to it. Standard output goes to *output_file* and standard error
    beside it with the suffix ``.err``. A log left by an earlier run is removed
    first, so it cannot be mistaken for this run's.

    The attempt's launch prefix (the parallel start, ``HTTK_WORKFLOW_LAUNCH``) is
    prepended by default; ``launch=False`` runs *argv* as given, and a command that
    already starts with a launcher such as ``srun`` or ``mpirun`` is refused with
    :class:`ValueError` when a prefix applies.

    :param argv: The LAMMPS command argument vector, without the input and log options.
    :param directory: Run LAMMPS in this directory.
    :param input_file: The input script name in *directory*.
    :param log_file: The log file name in *directory*.
    :param output_file: Save standard output under this name in *directory*.
    :param timeout: Stop the process after this many seconds when set.
    :param launch: Prepend the attempt's launch prefix when true, the default (``None``);
        ``False`` runs *argv* as given.
    :param termination_grace: Allow this many seconds for graceful termination.
    :param report_path: Write the report at this directory-relative path.
    :return: The classified run report.
    """

    root = Path(directory).resolve()
    (root / log_file).unlink(missing_ok=True)
    output = root / output_file
    # ponytail: no live monitor or remedy ladder; add them when a real campaign needs them.
    process = ProcessSupervisor().run(
        [*launch_command(argv, launch=launch is not False), "-in", input_file, "-log", log_file],
        timeout=timeout,
        cwd=root,
        termination_grace=termination_grace,
        stdout_path=output,
        stderr_path=output.with_suffix(".err"),
    )
    diagnostics = (*process.diagnostics, *diagnose_lammps(root, log=log_file, output=output_file))
    codes = {item.code for item in diagnostics}
    if process.timed_out:
        classification = "timeout"
    elif "lammps.error" in codes:
        classification = "crashed"
    elif process.returncode or "lammps.incomplete" in codes:
        classification = "process_failure"
    elif "lammps.minimization_not_converged" in codes:
        classification = "nonconverged"
    elif any(item.severity in {"error", "fatal"} for item in diagnostics):
        classification = "process_failure"
    else:
        classification = "completed"
    report = LammpsRunReport(process, classification, diagnostics, _parse(_read(root / log_file), _read(output)))
    report.write(root / report_path)
    return report
