"""The ``lammps`` code is registered through the ``codes`` registry tier, with its citation."""

import argparse
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

import httk.core  # noqa: F401  (importing httk.core runs registry discovery)
from httk.core.register import code_support, known_codes


def test_lammps_is_a_known_code_with_its_packaged_bash_api() -> None:
    assert "lammps" in known_codes()
    assert code_support("lammps").bash_api_path() == Path(str(files("httk.codes.lammps").joinpath("httk-lammps.sh")))


def test_the_bridge_mounts_the_lammps_commands() -> None:
    parser = argparse.ArgumentParser()
    code_support("lammps").resolve_bridge().add_commands(parser.add_subparsers(dest="command"))
    assert parser.parse_args(["lammps-thermo", "--column", "Temp"]).command == "lammps-thermo"


def test_the_lammps_credit_is_registered_on_import() -> None:
    script = """
from httk.core import credits
assert "Calculations with LAMMPS" not in credits.entries()
import httk.codes.lammps
assert len(credits.entries()["Calculations with LAMMPS"]) == 1
"""
    subprocess.run([sys.executable, "-c", script], check=True)
