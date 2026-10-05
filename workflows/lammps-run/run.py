#!/usr/bin/env python3
"""lammps.run: one LAMMPS run of one input script.

The single ``run`` step copies the ``script`` input (staged as
``files/in.lammps``) into the workdir, writes the optional ``structure`` input
(staged as ``files/POSCAR``) as the LAMMPS data file ``structure.data`` for a
script that says ``read_data structure.data``, runs LAMMPS under supervision,
and fails with the first diagnostic code (``lammps.error``, ...) when the run
is not clean.

Settings, resolved job parameter -> ``HTTK_*`` variable -> workspace setting:

* ``lammps.command``: the command that starts LAMMPS (default ``lmp``), the program
  only; the attempt's launch prefix supplies the parallel start.
"""

import shlex
import shutil

from httk.workflow import Attempt, Runner

from httk.codes.lammps import run_lammps, write_lammps_data

run = Runner("lammps.run")


@run.step(name="run")
def run_step(a: Attempt) -> None:
    """Prepare and run LAMMPS, then succeed or fail with what was diagnosed."""

    files = a.payload / "files"
    shutil.copyfile(files / "in.lammps", a.workdir / "in.lammps")
    if (files / "POSCAR").is_file():
        write_lammps_data(a.workdir / "structure.data", structure=files / "POSCAR")
    try:
        report = run_lammps(shlex.split(str(a.setting("lammps.command", "lmp"))), directory=a.workdir)
    except OSError as exception:
        a.fail("lammps.failed", f"could not start LAMMPS: {exception}")
        return
    if not report.ok:
        first = report.diagnostics[0] if report.diagnostics else None
        code = first.code if first else f"lammps.{report.classification}"
        a.fail(code, first.summary if first else f"LAMMPS {report.classification}")
        return
    a.state.merge({"thermo_last": report.result.thermo[-1].last if report.result.thermo else {}})
    a.succeed()


if __name__ == "__main__":
    raise SystemExit(run.main())
