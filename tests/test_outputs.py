"""``parse_lammps_log`` reads real captured LAMMPS 22 Jul 2025 (Update 6) logs."""

from conftest import DATA, LJ_LAST_ROW
from httk.codes.lammps import parse_lammps_log


def test_a_completed_md_run() -> None:
    result = parse_lammps_log(DATA / "lj.log")
    (table,) = result.thermo
    assert table.columns == ("Step", "Temp", "PotEng", "TotEng", "Press")
    assert len(table.rows) == 11 and table.rows[0] == (0.0, 1.0, -6.7733681, -5.2872569, -5.3989339)
    assert table.last == LJ_LAST_ROW
    assert (result.completed, result.errors, result.warnings, result.minimization_stop) == (True, (), (), None)
    assert result.minimization_converged is None


def test_two_runs_give_two_tables_past_interleaved_lines() -> None:
    result = parse_lammps_log(DATA / "two_runs.log")
    first, second = result.thermo
    # fix print lines interleave with the rows of both tables
    assert [row[0] for row in first.rows] == [0, 10, 20, 30, 40, 50]
    assert second.columns == ("Step", "Temp", "TotEng")
    assert second.last == {"Step": 80.0, "Temp": 0.54497405, "TotEng": -5.2940407}
    assert result.warnings == (
        "WARNING: New thermo_style command, previous thermo_modify settings will be lost (src/output.cpp:912)",
    )


def test_minimization_stopping_criteria() -> None:
    converged = parse_lammps_log(DATA / "min.log")
    assert (converged.minimization_stop, converged.minimization_converged) == ("energy tolerance", True)
    capped = parse_lammps_log(DATA / "min_maxiter.log")
    assert (capped.minimization_stop, capped.minimization_converged) == ("max iterations", False)


def test_an_error_stops_before_any_thermo_output() -> None:
    result = parse_lammps_log(DATA / "err.log")
    assert (result.thermo, result.completed) == ((), False)
    assert result.errors == ("ERROR: Unknown command: this_is_not_a_lammps_command foo bar (src/input.cpp:315)",)


def test_a_log_cut_short_keeps_the_rows_it_printed(tmp_path) -> None:
    text = (DATA / "lj.log").read_text(encoding="utf-8")
    (tmp_path / "log.lammps").write_text(text[: text.index("        60 ")], encoding="utf-8")
    result = parse_lammps_log(tmp_path / "log.lammps")
    assert result.completed is False and result.thermo[0].last["Step"] == 50.0
