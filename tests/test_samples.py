"""Synthetic dump/log joins pin run selection and unit conventions."""

from pathlib import Path

import pytest

pytest.importorskip("httk.atomistic")

from httk.core import definition_ids

from httk.codes.lammps import ThermoConversion, lammps_samples


def _sources(tmp_path: Path, *, units="metal", norm="no", steps=(0, 2, 4)):
    dump = tmp_path / "sample.lammpstrj"
    dump.write_text(
        "".join(
            f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n2\n"
            "ITEM: BOX BOUNDS pp pp pp\n0 4\n0 4\n0 4\n"
            "ITEM: ATOMS id type xu yu zu\n2 1 1 1 1\n1 1 0 0 0\n"
            for step in steps
        )
    )
    log = tmp_path / "log.lammps"
    log.write_text(
        f"units {units}\nthermo_style custom step pe press vol pxy v_test atoms\nthermo_modify norm {norm}\n"
        "run 4\nStep PotEng Press Volume Pxy v_test Atoms\n"
        "0 -2 10 64 2 3 2\n1 -3 10 64 2 3 2\n2 -4 10 64 2 3 2\n"
        "Loop time of 0.1 on 1 procs for 4 steps with 2 atoms\nTotal wall time: 0\n"
    )
    return dump, log


def test_intersection_converts_units_and_preserves_ids(tmp_path):
    dump, log = _sources(tmp_path)
    samples = tuple(
        lammps_samples(
            dump,
            log,
            species={1: "Ar"},
            units="metal",
            dimension=3,
            segment=0,
            table_index=0,
            columns=["PotEng", "Press", "Volume", "v_test"],
            observables=["atom_ids"],
            custom_columns={"v_test": ThermoConversion(definition_ids.TEMPERATURE, 2)},
        )
    )
    assert [sample.step for sample in samples] == [0, 2]
    first = dict(samples[0].thermo)
    assert first["potential_energy"] == -2
    assert first["volume"] == 64
    assert first["pressure"] == pytest.approx(0.001)
    assert first["temperature"] == 6
    assert [name for name, _ in samples[0].thermo] == ["potential_energy", "pressure", "volume", "temperature"]
    assert samples[0].observables == ((1, 2),)


def test_real_pressure_atm_to_gpa(tmp_path):
    dump, log = _sources(tmp_path, units="real")
    (row,) = tuple(
        lammps_samples(
            dump, log, species={1: "Ar"}, units="real", dimension=3, segment=0, table_index=0, columns=["Press"]
        )
    )[:1]
    assert row.thermo[0][1] == pytest.approx(10 * 101325 / 1e9)


def test_lj_units(tmp_path):
    dump, log = _sources(tmp_path, units="lj")
    log.write_text(
        log.read_text().replace(
            "Step PotEng Press Volume Pxy v_test Atoms", "Step PotEng Press Volume Temp v_test Atoms"
        )
    )
    (first, _) = tuple(
        lammps_samples(
            dump,
            log,
            species={1: "Ar"},
            units="lj",
            dimension=3,
            segment=0,
            table_index=0,
            columns=["PotEng", "Press", "Volume", "Temp"],
            trajectory_options={"lj_scales": (2, 1, 0.01)},
        )
    )
    got = dict(first.thermo)
    assert got["potential_energy"] == pytest.approx(-0.02)
    assert got["pressure"] == pytest.approx(10 * 0.01 / 8 * 160.2176634)
    assert got["volume"] == pytest.approx(64 * 8)
    assert got["temperature"] == pytest.approx(2 * 0.01 / 8.617333262145e-5)


def test_custom_stress_tensor_rejected(tmp_path):
    dump, log = _sources(tmp_path)
    with pytest.raises(ValueError, match="six"):
        tuple(
            lammps_samples(
                dump,
                log,
                species={1: "Ar"},
                units="metal",
                dimension=3,
                segment=0,
                table_index=0,
                columns=["v_test"],
                custom_columns={"v_test": ThermoConversion(definition_ids.STRESS_TENSOR, 1)},
            )
        )


def test_stress_tensor_is_one_tensile_voigt_value(tmp_path):
    dump, log = _sources(tmp_path)
    log.write_text(
        "units metal\nthermo_style custom step pxx pyy pzz pyz pxz pxy\nthermo_modify norm no\n"
        "run 4\nStep Pxx Pyy Pzz Pyz Pxz Pxy\n0 10 20 30 40 50 60\n2 10 20 30 40 50 60\n"
        "Loop time of 0.1 on 1 procs for 4 steps with 2 atoms\nTotal wall time: 0\n"
    )
    kwargs = {"species": {1: "Ar"}, "units": "metal", "dimension": 3, "segment": 0, "table_index": 0}
    with pytest.raises(ValueError, match="stress_basis"):
        tuple(lammps_samples(dump, log, columns=["Pxx", "Pyy", "Pzz", "Pyz", "Pxz", "Pxy"], **kwargs))
    (first, _) = tuple(
        lammps_samples(
            dump,
            log,
            columns=["Pxy", "Pxx", "Pyy", "Pzz", "Pyz", "Pxz"],
            stress_basis_matches_dump=True,
            **kwargs,
        )
    )
    assert len(first.thermo) == 1 and first.thermo[0][0] == "stress_tensor"
    assert first.thermo[0][1] == pytest.approx((-0.001, -0.002, -0.003, -0.004, -0.005, -0.006))


def test_strict_join_rejects_cadence_mismatch(tmp_path):
    dump, log = _sources(tmp_path)
    with pytest.raises(ValueError, match="step sequences"):
        tuple(
            lammps_samples(
                dump,
                log,
                species={1: "Ar"},
                units="metal",
                dimension=3,
                segment=0,
                table_index=0,
                columns=["PotEng"],
                join="strict",
            )
        )


def test_strict_join_exact_and_real_normalized_energy(tmp_path):
    dump, log = _sources(tmp_path, units="real", norm="yes", steps=(0, 1, 2))
    rows = tuple(
        lammps_samples(
            dump,
            log,
            species={1: "Ar"},
            units="real",
            dimension=3,
            segment=0,
            table_index=0,
            columns=["PotEng"],
            join="strict",
        )
    )
    assert rows[0].thermo[0][1] == pytest.approx(-4 * 4184 / (6.02214076e23 * 1.602176634e-19))


def test_normalized_global_count_cannot_be_inferred_from_dump_group(tmp_path):
    dump, log = _sources(tmp_path, norm="yes")
    text = log.read_text().replace("v_test Atoms", "v_test").replace("2 3 2\n", "2 3\n")
    log.write_text(text)
    with pytest.raises(ValueError, match="global atom count"):
        tuple(
            lammps_samples(
                dump, log, species={1: "Ar"}, units="metal", dimension=3, segment=0, table_index=0, columns=["PotEng"]
            )
        )
    rows = tuple(
        lammps_samples(
            dump,
            log,
            species={1: "Ar"},
            units="metal",
            dimension=3,
            segment=0,
            table_index=0,
            columns=["PotEng"],
            total_atom_count=10,
        )
    )
    assert rows[0].thermo[0][1] == -20


@pytest.mark.parametrize(
    "options,match",
    [
        ({"table_index": 4}, "table_index"),
        ({"columns": ["v_test"]}, "explicit ThermoConversion"),
        ({"columns": ["Pxy"]}, "all six"),
        ({"units": "real"}, "units"),
        ({"dimension": 2}, "dimension"),
        ({"join": "linear"}, "join"),
    ],
)
def test_invalid_join_metadata(tmp_path, options, match):
    dump, log = _sources(tmp_path)
    kwargs = {
        "species": {1: "Ar"},
        "units": "metal",
        "dimension": 3,
        "segment": 0,
        "table_index": 0,
        "columns": ["PotEng"],
    }
    kwargs.update(options)
    with pytest.raises(ValueError, match=match):
        tuple(lammps_samples(dump, log, **kwargs))
