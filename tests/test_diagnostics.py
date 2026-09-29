"""``diagnose_lammps`` maps finished runs to the stable ``lammps.*`` codes."""

import shutil
from pathlib import Path

from conftest import DATA
from httk.codes.lammps import diagnose_lammps


def _codes(directory: Path, log: str, output: str = "lammps.out") -> list[tuple[str, str]]:
    return [(item.code, item.severity) for item in diagnose_lammps(directory, log=log, output=output)]


def test_clean_runs_have_no_diagnostics() -> None:
    assert diagnose_lammps(DATA, log="lj.log", output="lj.out") == ()
    assert diagnose_lammps(DATA, log="two_runs.log") == ()  # a WARNING is not a diagnostic
    assert diagnose_lammps(DATA, log="min.log") == ()


def test_an_error_is_fatal_and_names_the_message() -> None:
    (diagnostic,) = diagnose_lammps(DATA, log="err.log", output="err.out")
    assert (diagnostic.code, diagnostic.severity) == ("lammps.error", "fatal")
    assert diagnostic.summary.startswith("ERROR: Unknown command: this_is_not_a_lammps_command")


def test_an_error_printed_only_to_the_screen_is_found(tmp_path: Path) -> None:
    # An unreadable input script: LAMMPS leaves an empty log and reports on stdout.
    shutil.copy(DATA / "noinput.out", tmp_path / "lammps.out")
    (tmp_path / "log.lammps").write_text("", encoding="utf-8")
    (diagnostic,) = diagnose_lammps(tmp_path)
    assert diagnostic.code == "lammps.error" and "Cannot open input script" in diagnostic.summary


def test_a_capped_minimization_is_not_converged() -> None:
    assert _codes(DATA, "min_maxiter.log") == [("lammps.minimization_not_converged", "error")]


def test_a_truncated_or_missing_log_is_incomplete(tmp_path: Path) -> None:
    text = (DATA / "lj.log").read_text(encoding="utf-8")
    (tmp_path / "log.lammps").write_text(text[: len(text) // 2], encoding="utf-8")
    assert _codes(tmp_path, "log.lammps") == [("lammps.incomplete", "error")]
    assert _codes(tmp_path, "absent.log") == [("lammps.incomplete", "error")]
