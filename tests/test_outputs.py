"""``parse_lammps_log`` reads real captured LAMMPS 22 Jul 2025 (Update 6) logs."""

import gzip

import pytest

from conftest import DATA, LJ_LAST_ROW
from httk.codes.lammps import KCAL_MOL_TO_EV, average_total_energy_ev, parse_lammps_log
from httk.codes.lammps.outputs import _parse


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
    assert [table.normalized for table in result.thermo] == [True, True]
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


def test_the_units_style_and_the_table_kinds_are_parsed() -> None:
    assert parse_lammps_log(DATA / "metal.log").units == "metal"
    assert parse_lammps_log(DATA / "metal.log").normalized is False
    assert parse_lammps_log(DATA / "real.log").normalized is False
    assert parse_lammps_log(DATA / "lj.log").normalized is True
    assert [table.kind for table in parse_lammps_log(DATA / "metal_min_run.log").thermo] == ["minimize", "run"]
    assert [table.kind for table in parse_lammps_log(DATA / "min.log").thermo] == ["minimize"]
    assert [table.kind for table in parse_lammps_log(DATA / "two_runs.log").thermo] == ["run", "run"]


def test_a_minimize_table_is_recognized_by_its_stats_block_without_echo() -> None:
    text = (DATA / "min.log").read_text(encoding="utf-8")
    result = _parse("\n".join(line for line in text.splitlines() if not line.startswith("minimize")))
    assert result.thermo[0].kind == "minimize"


def test_the_average_is_over_the_last_run_table_only() -> None:
    assert average_total_energy_ev(parse_lammps_log(DATA / "metal.log")) == -6.0
    assert average_total_energy_ev(parse_lammps_log(DATA / "metal_min_run.log")) == -6.0
    assert average_total_energy_ev(parse_lammps_log(DATA / "min.log")) is None
    assert average_total_energy_ev(parse_lammps_log(DATA / "metal_noetot.log")) is None


def test_real_units_convert_kcal_per_mol() -> None:
    assert KCAL_MOL_TO_EV == pytest.approx(0.0433641043, rel=1e-8)
    assert average_total_energy_ev(parse_lammps_log(DATA / "real.log")) == pytest.approx(-6.0 * KCAL_MOL_TO_EV)


def test_unsupported_units_are_refused() -> None:
    text = (
        (DATA / "lj.log")
        .read_text(encoding="utf-8")
        .replace(
            "thermo_style custom step temp pe etotal press",
            "thermo_style custom step temp pe etotal press\nthermo_modify norm no",
        )
    )
    with pytest.raises(ValueError, match="'lj'"):
        average_total_energy_ev(_parse(text))


def test_a_compressed_log_is_parsed(tmp_path) -> None:
    path = tmp_path / "log.lammps.gz"
    path.write_bytes(gzip.compress((DATA / "metal.log").read_bytes()))
    assert average_total_energy_ev(parse_lammps_log(path)) == -6.0


def test_a_trailing_run_0_does_not_replace_the_md_table() -> None:
    text = (DATA / "metal.log").read_text(encoding="utf-8")
    first = text.split("Loop time of", 1)[0]
    table = first[first.index("   Step") :]
    zero = "run 0\n" + table.split("\n", 2)[0] + "\n" + table.split("\n", 2)[1].replace("-1.0", "-99.0") + "\n"
    zero += "Loop time of 0.1 on 1 procs for 0 steps with 108 atoms\n"
    result = _parse(text + "\n" + zero)
    assert [table.steps for table in result.thermo] == [100, 0]
    assert average_total_energy_ev(result) == -6.0


def test_thermo_modify_norm_is_detected_and_refused() -> None:
    result = parse_lammps_log(DATA / "two_runs.log")
    assert result.normalized
    metal = _parse(
        (DATA / "metal.log").read_text(encoding="utf-8").replace("run         100", "thermo_modify norm yes\nrun 100")
    )
    with pytest.raises(ValueError, match="per-atom"):
        average_total_energy_ev(metal)


def test_thermo_style_resets_norm_to_the_units_default() -> None:
    text = (DATA / "metal.log").read_text(encoding="utf-8")
    text = text.replace(
        "thermo_style custom step temp pe etotal press",
        "thermo_modify norm yes\nthermo_style custom step temp pe etotal press",
    )
    result = _parse(text)
    assert result.thermo[0].normalized is False


def test_unknown_norm_state_is_preserved_and_refused() -> None:
    text = (DATA / "metal.log").read_text(encoding="utf-8")
    without_units = text.replace("units       metal\n", "")
    result = _parse(without_units)
    assert result.thermo[0].normalized is None
    with pytest.raises(ValueError, match="normalization state is unknown"):
        average_total_energy_ev(result)
    after_clear = _parse(text.replace("units       metal", "units metal\nclear", 1))
    assert after_clear.thermo[0].normalized is None


def test_each_thermo_table_keeps_its_own_norm_state() -> None:
    text = (DATA / "two_runs.log").read_text(encoding="utf-8").replace("units       lj", "units metal", 1)
    text = text.replace(
        "thermo_style custom step temp etotal\n", "thermo_style custom step temp etotal\nthermo_modify norm yes\n"
    )
    result = _parse(text)
    assert [table.normalized for table in result.thermo] == [False, True]
    assert result.normalized is True


def test_invalid_norm_value_does_not_become_false() -> None:
    text = (
        (DATA / "metal.log").read_text(encoding="utf-8").replace("run         100", "thermo_modify norm maybe\nrun 100")
    )
    result = _parse(text)
    assert result.thermo[0].normalized is None


def _metal(*, before: str = "", replace_units: str = "units       metal") -> str:
    text = (DATA / "metal.log").read_text(encoding="utf-8")
    return text.replace("units       metal", before + replace_units, 1)


def test_the_units_in_force_at_the_table_are_used() -> None:
    cleared = _parse(_metal(before="units real\nclear\n"))
    assert average_total_energy_ev(cleared) == -6.0
    twice = _parse(_metal(before="units real\n"))
    assert average_total_energy_ev(twice) == -6.0
    assert twice.thermo[0].units == "metal"
    real_last = _parse(_metal(before="units metal\n", replace_units="units real"))
    assert average_total_energy_ev(real_last) == pytest.approx(-6.0 * KCAL_MOL_TO_EV)


def test_norm_yes_then_no_is_accepted() -> None:
    text = (DATA / "metal.log").read_text(encoding="utf-8")
    both = _parse(text.replace("run         100", "thermo_modify norm yes\nthermo_modify lost no norm no\nrun 100"))
    assert not both.thermo[0].normalized
    assert average_total_energy_ev(both) == -6.0
