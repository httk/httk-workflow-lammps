"""Shared test configuration: isolate the httk configuration of every test."""

import os
import shlex
import shutil
from pathlib import Path

import pytest

# Keep each BLAS/OpenMP runtime of the many short-lived runner processes to one
# thread; child runners (and LAMMPS) inherit this.
for _thread_limit in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_thread_limit] = "1"


@pytest.fixture(autouse=True)
def _isolated_httk_config(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    """Give every test its own httk config and data home, so no workspace registry leaks between tests."""

    monkeypatch.setenv("HTTK_CONFIG_HOME", str(tmp_path_factory.mktemp("httk-config")))
    monkeypatch.setenv("HTTK_DATA_HOME", str(tmp_path_factory.mktemp("httk-store")))
    # A developer's launch prefix must not leak into the tests.
    monkeypatch.delenv("HTTK_WORKFLOW_LAUNCH", raising=False)


DATA = Path(__file__).resolve().parent / "data"
REPO_ROOT = Path(__file__).resolve().parent.parent

#: The final thermo row of the captured tests/data/lj.log (Step Temp PotEng TotEng Press).
LJ_LAST_ROW = {"Step": 100.0, "Temp": 0.57389677, "PotEng": -6.1447617, "TotEng": -5.2918874, "Press": -2.0049305}

# Conventional fcc argon (a = 5.26 Å), for a script with an LJ pair style.
ARGON_POSCAR = """fcc argon
5.26
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Ar
4
Direct
0.0 0.0 0.0
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
"""

# Wurtzite-like hexagonal cell: a = 3.2 Å, c = 5.2 Å, gamma = 120 degrees.
HEXAGONAL_POSCAR = """hexagonal
1.0
3.2 0.0 0.0
-1.6 2.7712812921102 0.0
0.0 0.0 5.2
Zn O
2 2
Direct
0.3333333333 0.6666666667 0.0
0.6666666667 0.3333333333 0.5
0.3333333333 0.6666666667 0.38
0.6666666667 0.3333333333 0.88
"""

# Reads ARGON_POSCAR written as structure.data and evaluates it once.
ARGON_SCRIPT = """units metal
atom_style atomic
boundary p p p
read_data structure.data
mass 1 39.948
pair_style lj/cut 8.0
pair_coeff 1 1 0.0104 3.40
thermo_style custom step pe press
run 0
"""


def lammps_command() -> list[str] | None:
    """The real LAMMPS command: ``HTTK_TEST_LAMMPS_COMMAND``, else ``lmp`` on PATH, else ``None``."""

    command = os.environ.get("HTTK_TEST_LAMMPS_COMMAND") or shutil.which("lmp")
    return shlex.split(command) if command else None


requires_lammps = pytest.mark.skipif(
    lammps_command() is None, reason="needs a real LAMMPS: set HTTK_TEST_LAMMPS_COMMAND or put lmp on PATH"
)
