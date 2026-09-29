"""``run_lammps`` classifies a supervised run and writes its report.

A stand-in ``lmp`` replays a captured log, so the classification is exercised
without a real LAMMPS.
"""

import json
import sys
from pathlib import Path

import pytest

from conftest import DATA, LJ_LAST_ROW
from httk.codes.lammps import run_lammps

# Copies the named captured log to the -log path run_lammps appends and exits with a code.
_REPLAY = "import shutil, sys; shutil.copy(sys.argv[1], sys.argv[-1]); sys.exit(int(sys.argv[2]))"


def _replay(log: str, code: int = 0) -> list[str]:
    return [sys.executable, "-c", _REPLAY, str(DATA / log), str(code)]


@pytest.mark.parametrize(
    ("argv", "classification"),
    [
        (_replay("lj.log"), "completed"),
        (_replay("min_maxiter.log"), "nonconverged"),
        (_replay("err.log", 1), "crashed"),
        (_replay("lj.log", 3), "process_failure"),
    ],
)
def test_the_run_is_classified_and_reported(tmp_path: Path, argv: list[str], classification: str) -> None:
    report = run_lammps(argv, directory=tmp_path)
    assert report.classification == classification
    assert report.ok == (classification == "completed")
    assert report.process.argv[-4:] == ("-in", "in.lammps", "-log", "log.lammps")
    saved = json.loads((tmp_path / "lammps-run-report.json").read_text(encoding="utf-8"))
    assert saved["format"] == "httk-lammps-run-report" and saved["classification"] == classification
    if classification == "completed":
        assert saved["result"]["thermo"][-1]["rows"][-1] == list(LJ_LAST_ROW.values())
        assert report.diagnostics == ()


def test_a_stale_log_is_removed_before_the_run(tmp_path: Path) -> None:
    (tmp_path / "log.lammps").write_text((DATA / "lj.log").read_text(encoding="utf-8"), encoding="utf-8")
    report = run_lammps([sys.executable, "-c", "pass"], directory=tmp_path)
    assert report.classification == "process_failure"
    assert [item.code for item in report.diagnostics] == ["lammps.incomplete"]


def test_a_clean_run_after_an_error_in_the_same_directory_is_completed(tmp_path: Path) -> None:
    # A stand-in lmp that prints a captured stdout and copies a captured log to the -log path.
    replay = "import shutil, sys; print(open(sys.argv[1]).read()); shutil.copy(sys.argv[2], sys.argv[-1])"
    lmp = [sys.executable, "-c", replay]
    failed = run_lammps([*lmp, str(DATA / "err.out"), str(DATA / "err.log")], directory=tmp_path)
    assert failed.classification == "crashed"
    clean = run_lammps([*lmp, str(DATA / "lj.out"), str(DATA / "lj.log")], directory=tmp_path)
    assert clean.classification == "completed", clean.result.errors
