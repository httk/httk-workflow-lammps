"""``write_lammps_data`` writes a data file the real LAMMPS reads."""

import logging
import subprocess
from pathlib import Path

import pytest

from conftest import ARGON_POSCAR, ARGON_SCRIPT, HEXAGONAL_POSCAR, lammps_command, requires_lammps
from httk.codes.lammps import parse_lammps_log, write_lammps_data

pytest.importorskip("httk.atomistic")

# HEXAGONAL_POSCAR turned 90 degrees about z: the same crystal in another frame.
ROTATED_HEXAGONAL_POSCAR = HEXAGONAL_POSCAR.replace(
    "3.2 0.0 0.0\n-1.6 2.7712812921102 0.0", "0.0 3.2 0.0\n-2.7712812921102 -1.6 0.0"
)


def _write(tmp_path: Path, poscar: str, **options: object) -> list[str]:
    (tmp_path / "POSCAR").write_text(poscar, encoding="utf-8")
    write_lammps_data(tmp_path / "structure.data", structure=tmp_path / "POSCAR", **options)  # type: ignore[arg-type]
    return (tmp_path / "structure.data").read_text(encoding="utf-8").splitlines()


def _numbers(lines: list[str], start: str, count: int) -> list[float]:
    index = lines.index(start) + 2
    return [float(x) for line in lines[index : index + count] for x in line.split("#")[0].split()]


def test_a_cubic_cell_is_an_orthogonal_box(tmp_path: Path) -> None:
    lines = _write(tmp_path, ARGON_POSCAR)
    assert lines[2:8] == ["4 atoms", "1 atom types", "", "0.0 5.26 xlo xhi", "0.0 5.26 ylo yhi", "0.0 5.26 zlo zhi"]
    assert not any(line.endswith("xy xz yz") for line in lines) and "Masses" not in lines
    assert "Atoms # atomic" in lines
    assert _numbers(lines, "Atoms # atomic", 4) == pytest.approx(
        [1, 1, 0, 0, 0, 2, 1, 0, 2.63, 2.63, 3, 1, 2.63, 0, 2.63, 4, 1, 2.63, 2.63, 0]
    )


@pytest.mark.parametrize("poscar", [HEXAGONAL_POSCAR, ROTATED_HEXAGONAL_POSCAR])
def test_a_hexagonal_cell_is_a_restricted_triclinic_box(tmp_path: Path, poscar: str) -> None:
    lines = _write(tmp_path, poscar, masses={"Zn": 65.38, "O": 15.999})
    assert lines[0].endswith("types: Zn O") and lines[3] == "2 atom types"
    box = [float(x) for line in lines[5:9] for x in line.split()[:3] if x[0].isdigit() or x[0] == "-"]
    assert box == pytest.approx([0, 3.2, 0, 2.7712812921102, 0, 5.2, -1.6, 0, 0], abs=1e-9)
    assert lines[8].endswith("xy xz yz")
    assert _numbers(lines, "Masses", 2) == [1, 65.38, 2, 15.999]
    positions = _numbers(lines, "Atoms # atomic", 4)
    assert positions[:5] == pytest.approx([1, 1, 0, 1.8475208614, 0], abs=1e-6)
    assert positions[10:15] == pytest.approx([3, 2, 0, 1.8475208614, 1.976], abs=1e-6)


def test_a_poscar_is_written_as_the_floats_it_writes(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    # As an exact structure, -1.92 would become -23/12, the simplest fraction within its digits.
    poscar = HEXAGONAL_POSCAR.replace("3.2 0.0 0.0\n-1.6 2.7712812921102", "3.84 0.0 0.0\n-1.92 3.33")
    with caplog.at_level(logging.WARNING):
        lines = _write(tmp_path, poscar)
    assert caplog.records == []
    assert lines[5:9] == ["0.0 3.84 xlo xhi", "0.0 3.33 ylo yhi", "0.0 5.2 zlo zhi", "-1.92 0.0 0.0 xy xz yz"]
    # 0.3333333333 a + 0.6666666667 b as written, not a + 2b over 3
    assert lines[-4] == "1 1 -1.92e-10 2.220000000111 0.0"


def test_a_structure_object_is_written_from_its_exact_values(tmp_path: Path) -> None:
    import httk.core

    (tmp_path / "POSCAR").write_text(ARGON_POSCAR, encoding="utf-8")
    structure = httk.core.load(str(tmp_path / "POSCAR"), precision=1e-6)
    write_lammps_data(tmp_path / "structure.data", structure=structure)
    lines = (tmp_path / "structure.data").read_text(encoding="utf-8").splitlines()
    assert lines[5] == "0.0 5.26 xlo xhi" and lines[-1] == "4 1 2.63 2.63 0.0"


def test_invalid_options_are_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="no mass given for species O"):
        _write(tmp_path, HEXAGONAL_POSCAR, masses={"Zn": 65.38})
    with pytest.raises(ValueError, match="atom_style"):
        _write(tmp_path, ARGON_POSCAR, atom_style="charge")
    left_handed = ARGON_POSCAR.replace("1.0 0.0 0.0\n0.0 1.0 0.0", "0.0 1.0 0.0\n1.0 0.0 0.0")
    with pytest.raises(ValueError, match="left-handed"):
        _write(tmp_path, left_handed)


@requires_lammps
def test_the_real_lammps_reads_the_data_file(tmp_path: Path) -> None:
    command = lammps_command()
    assert command is not None
    _write(tmp_path, HEXAGONAL_POSCAR, masses={"Zn": 65.38, "O": 15.999})
    script = ARGON_SCRIPT.replace("mass 1 39.948\n", "").replace("pair_coeff 1 1", "pair_coeff * *")
    (tmp_path / "in.lammps").write_text(script, encoding="utf-8")
    subprocess.run(
        [*command, "-in", "in.lammps", "-log", "log.lammps"],
        cwd=tmp_path,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=True,
        timeout=120,
    )
    log = (tmp_path / "log.lammps").read_text(encoding="utf-8")
    assert "  triclinic box = (0 0 0) to (3.2 2.7712813 5.2) with tilt (-1.6 0 0)" in log
    assert "4 atoms" in log
    result = parse_lammps_log(tmp_path / "log.lammps")
    assert result.completed and result.errors == () and result.thermo[0].columns == ("Step", "PotEng", "Press")
