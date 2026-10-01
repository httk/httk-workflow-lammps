"""The ``lammps.calculation`` collector recognizes and collects free-standing runs via ``collect_tree``."""

import bz2
import shutil
from pathlib import Path

import pytest
from httk.workflow import claims, collect_tree
from httk.workflow.calculations import content_digest

from conftest import DATA
from httk.codes.lammps import KCAL_MOL_TO_EV
from httk.codes.lammps.collect import find_outputs, read_average_total_energy

SCRIPT = "units metal\nrun 100\n"


def _run(directory: Path, log: str = "metal.log", name: str = "log.lammps", script: str | None = "in.lammps") -> Path:
    directory.mkdir(parents=True)
    shutil.copy(DATA / log, directory / name)
    if script:
        (directory / script).write_text(SCRIPT, encoding="utf-8")
    return directory


def test_a_metal_run_is_collected(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    (item,) = collect_tree(tmp_path)
    assert item.missing_collector is None
    assert item.run.source_id == f"lammps.calculation:{content_digest(directory, ['in.lammps'])}"
    assert item.outputs["average_total_energy"].value == -6.0  # type: ignore[attr-defined]


def test_a_real_run_converts_kcal_per_mol(tmp_path: Path) -> None:
    _run(tmp_path / "md", "real.log", "run.log")
    (item,) = collect_tree(tmp_path)
    assert item.outputs["average_total_energy"].value == pytest.approx(-6.0 * KCAL_MOL_TO_EV)  # type: ignore[attr-defined]


def test_the_marker_needs_a_dot_lammps_file(tmp_path: Path) -> None:
    _run(tmp_path / "md", "metal.log", "run.log", "md.in")
    assert list(claims(tmp_path)) == []


def test_lj_units_are_unclaimed(tmp_path: Path) -> None:
    _run(tmp_path / "md", "lj.log")
    (outcome,) = claims(tmp_path)
    assert (outcome.kind, outcome.reason) == ("unclaimed", "LAMMPS units lj cannot be converted to eV")
    assert list(collect_tree(tmp_path)) == []


def test_only_the_last_run_table_is_averaged(tmp_path: Path) -> None:
    _run(tmp_path / "md", "metal_min_run.log")
    (item,) = collect_tree(tmp_path)
    assert item.outputs["average_total_energy"].value == -6.0  # type: ignore[attr-defined]


def test_a_log_without_toteng_is_claimed_and_degraded(tmp_path: Path) -> None:
    _run(tmp_path / "md", "metal_noetot.log")
    (item,) = collect_tree(tmp_path)
    assert item.missing_collector is not None and "print etotal in thermo_style" in item.missing_collector


def test_an_incomplete_run_is_degraded(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    log = directory / "log.lammps"
    log.write_text(log.read_text(encoding="utf-8").replace("Total wall time", "Total wall"), encoding="utf-8")
    (item,) = collect_tree(tmp_path)
    assert item.missing_collector is not None and "not a completed LAMMPS run" in item.missing_collector


def test_a_missing_or_ambiguous_script_is_unclaimed(tmp_path: Path) -> None:
    directory = _run(tmp_path / "a", script=None)
    (outcome,) = claims(tmp_path)
    assert (
        outcome.kind == "unclaimed"
        and outcome.reason
        and "no in.* or *.in input script beside log.lammps" in outcome.reason
    )
    (directory / "in.a").write_text(SCRIPT, encoding="utf-8")
    (directory / "b.in").write_text(SCRIPT, encoding="utf-8")
    (outcome,) = claims(tmp_path)
    assert outcome.reason == "several LAMMPS input scripts: b.in, in.a"


def test_several_logs_are_unclaimed(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    (directory / "log.lammps").rename(directory / "x.lammps")
    shutil.copy(DATA / "metal.log", directory / "y.lammps")
    (outcome,) = claims(tmp_path)
    assert outcome.reason == "several LAMMPS logs: x.lammps, y.lammps"


def test_a_captured_screen_copy_beside_the_log_is_ignored(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    shutil.copy(DATA / "lj.out", directory / "slurm-1.out")
    (item,) = collect_tree(tmp_path)
    assert item.outputs["average_total_energy"].value == -6.0  # type: ignore[attr-defined]
    assert [path.name for path in find_outputs(directory)] == ["log.lammps"]


def test_a_log_without_echoed_input_is_unclaimed(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    log = directory / "log.lammps"
    log.write_text(log.read_text(encoding="utf-8").replace("units       metal\n", ""), encoding="utf-8")
    (outcome,) = claims(tmp_path)
    assert (outcome.kind, outcome.reason) == ("unclaimed", "LAMMPS units unknown cannot be converted to eV")


def test_norm_yes_is_unclaimed_and_refused_by_the_reader(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    log = directory / "log.lammps"
    log.write_text(log.read_text(encoding="utf-8").replace("run         100", "thermo_modify norm yes\nrun 100"))
    (outcome,) = claims(tmp_path)
    assert outcome.kind == "unclaimed" and outcome.reason and "per-atom" in outcome.reason
    with pytest.raises(ValueError, match="per-atom"):
        read_average_total_energy(log)


def test_files_read_by_the_script_are_part_of_the_identity(tmp_path: Path) -> None:
    for name, data in (("a", "box 1\n"), ("b", "box 2\n")):
        directory = _run(tmp_path / name)
        (directory / "in.lammps").write_text(SCRIPT + "read_data data.lmp\n", encoding="utf-8")
        (directory / "data.lmp").write_text(data, encoding="utf-8")
    items = list(collect_tree(tmp_path))
    assert len(items) == 2
    assert len({item.run.source_id for item in items}) == 2


def test_a_compressed_log_is_collected(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    (directory / "log.lammps.bz2").write_bytes(bz2.compress((directory / "log.lammps").read_bytes()))
    (directory / "log.lammps").unlink()
    (item,) = collect_tree(tmp_path)
    assert item.outputs["average_total_energy"].value == -6.0  # type: ignore[attr-defined]


def test_a_scheduler_log_is_no_candidate(tmp_path: Path) -> None:
    (tmp_path / "slurm-1.out").write_text("Job started\n" * 5, encoding="utf-8")
    (tmp_path / "run.log").write_text("not LAMMPS\n", encoding="utf-8")
    (tmp_path / "in.lammps").write_text(SCRIPT, encoding="utf-8")
    assert list(claims(tmp_path)) == []
    assert find_outputs(tmp_path) == ()


def test_a_different_potential_file_changes_the_identity(tmp_path: Path) -> None:
    for name, potential in (("a", "Fe 1\n"), ("b", "Fe 2\n")):
        directory = _run(tmp_path / name)
        (directory / "in.lammps").write_text(SCRIPT + "pair_coeff * * Fe.eam.alloy Fe\n", encoding="utf-8")
        (directory / "Fe.eam.alloy").write_text(potential, encoding="utf-8")
    items = list(collect_tree(tmp_path))
    assert len({item.run.source_id for item in items}) == 2


def test_the_units_after_a_clear_decide_the_claim(tmp_path: Path) -> None:
    directory = _run(tmp_path / "md")
    log = directory / "log.lammps"
    log.write_text(log.read_text(encoding="utf-8").replace("units       metal", "units real\nclear\nunits metal"))
    (item,) = collect_tree(tmp_path)
    assert item.outputs["average_total_energy"].value == -6.0  # type: ignore[attr-defined]
