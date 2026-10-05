"""The ``lammps-*`` bridge commands and the Bash API that forwards to them.

Code verbs need no attempt, so they run straight through the shell bridge.
"""

import json
import os
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

import pytest

from conftest import DATA, HEXAGONAL_POSCAR


def _bridge(cwd: Path, *arguments: str) -> "subprocess.CompletedProcess[str]":
    environment = {name: value for name, value in os.environ.items() if not name.startswith("HTTK_WORKFLOW_")}
    environment["PYTHONPATH"] = str(Path(__file__).parents[1] / "src")
    return subprocess.run(
        [sys.executable, "-m", "httk.workflow._shell_bridge", *arguments],
        cwd=cwd,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )


def test_thermo_answers_and_absences(tmp_path: Path) -> None:
    last = _bridge(tmp_path, "lammps-thermo", "--log", str(DATA / "lj.log"), "--column", "TotEng")
    assert (last.returncode, last.stdout) == (0, "-5.2918874\n")
    first = _bridge(tmp_path, "lammps-thermo", "--log", str(DATA / "two_runs.log"), "--column", "Step", "--table", "0")
    assert (first.returncode, first.stdout) == (0, "50\n")
    for extra in (["--column", "Press", "--table", "1"], ["--column", "Nope"], ["--column", "Step", "--table", "5"]):
        absent = _bridge(tmp_path, "lammps-thermo", "--log", str(DATA / "two_runs.log"), *extra)
        assert (absent.returncode, absent.stdout) == (1, "")
    no_tables = _bridge(tmp_path, "lammps-thermo", "--log", str(DATA / "err.log"), "--column", "Step")
    assert (no_tables.returncode, no_tables.stdout) == (1, "")
    refused = _bridge(tmp_path, "lammps-thermo", "--log", str(tmp_path / "missing.log"), "--column", "Step")
    assert refused.returncode == 2


def test_diagnose_prints_codes_and_json(tmp_path: Path) -> None:
    clean = _bridge(tmp_path, "lammps-diagnose", "--log", str(DATA / "lj.log"), "--output", "lj.out")
    assert (clean.returncode, clean.stdout) == (0, "")
    error = _bridge(tmp_path, "lammps-diagnose", "--log", str(DATA / "err.log"), "--output", "err.out", "--json")
    assert error.returncode == 20
    assert [item["code"] for item in json.loads(error.stdout)] == ["lammps.error"]


def test_write_data_then_run_with_a_replayed_lmp(tmp_path: Path) -> None:
    pytest.importorskip("httk.atomistic")
    (tmp_path / "POSCAR").write_text(HEXAGONAL_POSCAR, encoding="utf-8")
    options = {"structure": "POSCAR", "masses": {"Zn": 65.38, "O": 15.999}}
    (tmp_path / "options.json").write_text(json.dumps(options), encoding="utf-8")
    assert _bridge(tmp_path, "lammps-write-data", "--options", "options.json").returncode == 0
    assert "2 atom types" in (tmp_path / "structure.data").read_text(encoding="utf-8")

    replay = "import shutil, sys; shutil.copy(sys.argv[1], sys.argv[-1])"
    ran = _bridge(tmp_path, "lammps-run", "--", sys.executable, "-c", replay, str(DATA / "min_maxiter.log"))
    assert (ran.returncode, ran.stdout) == (21, "lammps-run-report.json\n")
    report = json.loads((tmp_path / "lammps-run-report.json").read_text(encoding="utf-8"))
    assert report["classification"] == "nonconverged"


def test_the_bash_api_forwards_to_the_bridge(tmp_path: Path) -> None:
    workflow_api = files("httk.workflow").joinpath("languages", "bash", "httk-workflow.sh")
    lammps_api = files("httk.codes.lammps").joinpath("httk-lammps.sh")
    script = (
        f'source "{workflow_api}"; source "{lammps_api}"; '
        f'httk_lammps_thermo --log "{DATA / "lj.log"}" --column Temp; echo "v$HTTK_LAMMPS_BASH_API_VERSION"'
    )
    environment = {name: value for name, value in os.environ.items() if not name.startswith("HTTK_WORKFLOW_")}
    environment["PYTHONPATH"] = str(Path(__file__).parents[1] / "src")
    environment["HTTK_WORKFLOW_PYTHON"] = sys.executable
    result = subprocess.run(
        ["bash", "-c", script], cwd=tmp_path, env=environment, text=True, capture_output=True, check=False
    )
    assert (result.returncode, result.stdout) == (0, "0.57389677\nv1\n"), result.stderr
    unguarded = subprocess.run(
        ["bash", "-c", f'source "{lammps_api}"; httk_lammps_thermo --column Temp'],
        text=True,
        capture_output=True,
        check=False,
    )
    assert unguarded.returncode == 2 and "source HTTK_WORKFLOW_BASH_API" in unguarded.stderr


def test_no_launch_runs_the_command_as_given(tmp_path: Path) -> None:
    # _bridge strips HTTK_WORKFLOW_*; call the bridge directly with a launch prefix set.
    environment = {
        **os.environ,
        "PYTHONPATH": str(Path(__file__).parents[1] / "src"),
        "HTTK_WORKFLOW_LAUNCH": "env A=b",
    }
    command = [sys.executable, "-m", "httk.workflow._shell_bridge", "lammps-run"]
    program = ["--", sys.executable, "-c", "pass"]

    def run(*flags: str) -> "subprocess.CompletedProcess[str]":
        return subprocess.run(
            [*command, *flags, *program], cwd=tmp_path, env=environment, text=True, capture_output=True, check=False
        )

    assert run("--no-launch").returncode == 22
    report = json.loads((tmp_path / "lammps-run-report.json").read_text(encoding="utf-8"))
    assert report["process"]["argv"][0] == sys.executable
    assert run().returncode == 22
    report = json.loads((tmp_path / "lammps-run-report.json").read_text(encoding="utf-8"))
    assert report["process"]["argv"][:2] == ["env", "A=b"]
