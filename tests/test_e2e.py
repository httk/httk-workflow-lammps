"""End to end with the real LAMMPS: install, run, and collect the ``lammps.run`` workflow.

Runs only with a real LAMMPS (``HTTK_TEST_LAMMPS_COMMAND`` or ``lmp`` on PATH).
The workflow is installed as the ``httk_plugin.toml`` plugin of this repository,
as ``httk plugin install`` does, into the test's isolated data home.
"""

import json
import shlex
import shutil
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

from conftest import ARGON_POSCAR, ARGON_SCRIPT, DATA, LJ_LAST_ROW, REPO_ROOT, lammps_command, requires_lammps
from httk.codes.lammps import parse_lammps_log

pytestmark = [requires_lammps, pytest.mark.slow]


@pytest.fixture
def installed_plugin(tmp_path_factory: pytest.TempPathFactory) -> Iterator[None]:
    from httk.core.plugins import install_plugin
    from httk.workflow.packages import _reset_plugin_workflow_cache

    source = tmp_path_factory.mktemp("plugin-source") / "httk-workflow-lammps"
    shutil.copytree(REPO_ROOT / "workflows", source / "workflows", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy2(REPO_ROOT / "httk_plugin.toml", source)
    install_plugin(source)
    _reset_plugin_workflow_cache()
    yield
    _reset_plugin_workflow_cache()


def _run_until_idle(workspace: Any, job: Any) -> Path:
    from httk.workflow import TaskManager

    with TaskManager(workspace, heartbeat_interval=0.01) as manager:
        manager.run_until_idle(timeout=600.0)
    marker = workspace.find_marker_by_id(job.job_id)
    assert marker is not None
    assert marker.kind == "succeeded", workspace.read_state(marker).get("failure")
    (log,) = workspace.root.rglob("log.lammps")
    return log


def test_lammps_run_melts_lj_and_collects_a_run_only_record(
    tmp_path: Path, installed_plugin: None, capsys: pytest.CaptureFixture[str]
) -> None:
    store = pytest.importorskip("httk.store")
    from httk.core import DataRecord, Run
    from httk.core.cli import CLIContext
    from httk.workflow import Workspace
    from httk.workflow.registry import register_workspace
    from httk.workflow.scaffold import new_job
    from httk.workflow.workflow_cli import command

    workspace = Workspace.initialize(tmp_path / "workspace")
    workspace.set_setting("lammps.command", shlex.join(lammps_command() or ()))
    job = new_job(workspace, "lammps.run", inputs={"script": DATA / "lj.in"})
    last = parse_lammps_log(_run_until_idle(workspace, job)).thermo[-1].last
    assert last.keys() == LJ_LAST_ROW.keys()
    assert list(last.values()) == pytest.approx(list(LJ_LAST_ROW.values()), abs=1e-6)

    register_workspace("lammps", str(workspace.root))
    database = tmp_path / "results.sqlite"
    arguments = ["collect", "--workspace", "lammps", "--into", str(database), "--id-base", "httk.test"]
    assert command([*arguments, "--no-id-ledger"], CLIContext("httk", tmp_path)) == 0
    report = json.loads(capsys.readouterr().out.splitlines()[0])
    assert report["run_only"] is True and report["outputs"] in ({}, []) and report["stored"]["run"]

    with store.Backend.sqlite(database) as backend:
        searcher = store.SqlStore(backend).searcher()
        runs = list(searcher.results(run=searcher.variable(Run)))
        searcher = store.SqlStore(backend).searcher()
        records = list(searcher.results(record=searcher.variable(DataRecord)))
    assert len(runs) == 1 and records == []


def test_lammps_run_reads_a_structure_written_as_a_data_file(tmp_path: Path, installed_plugin: None) -> None:
    pytest.importorskip("httk.atomistic")
    from httk.workflow import Workspace
    from httk.workflow.scaffold import new_job

    workspace = Workspace.initialize(tmp_path / "workspace")
    workspace.set_setting("lammps.command", shlex.join(lammps_command() or ()))
    (tmp_path / "in.lammps").write_text(ARGON_SCRIPT, encoding="utf-8")
    (tmp_path / "POSCAR").write_text(ARGON_POSCAR, encoding="utf-8")
    job = new_job(workspace, "lammps.run", inputs={"script": tmp_path / "in.lammps", "structure": tmp_path / "POSCAR"})
    log = _run_until_idle(workspace, job)
    assert "  4 atoms" in log.read_text(encoding="utf-8")
    (table,) = parse_lammps_log(log).thermo
    # fcc argon near its LJ equilibrium: bound, about -0.08 eV per atom for 4 atoms (pe is not normalized)
    assert table.columns == ("Step", "PotEng", "Press") and -0.4 < table.last["PotEng"] < -0.2
