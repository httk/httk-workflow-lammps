"""Synthetic dump/log joins pin run selection and unit conventions."""

from pathlib import Path

import pytest

pytest.importorskip("httk.atomistic")

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
            columns=["PotEng", "Press", "Volume", "Pxy", "v_test"],
            observables=["atom_ids"],
            stress_basis_matches_dump=True,
            custom_columns={"v_test": ThermoConversion("custom", "K", 2)},
        )
    )
    assert [sample.step for sample in samples] == [0, 2]
    first = {name: value for name, value, _ in samples[0].thermo}
    assert first["potential_energy"] == -2
    assert first["volume"] == 64
    assert first["pressure"] == pytest.approx(10e5 / 1.602176634e11)
    assert first["stress_xy"] == pytest.approx(-2e5 / 1.602176634e11)
    assert first["custom"] == 6
    assert samples[0].observables == ((1, 2),)


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
        ({"columns": ["Pxy"]}, "stress_basis"),
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
