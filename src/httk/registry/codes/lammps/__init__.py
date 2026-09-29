"""Register the LAMMPS code support implemented by :mod:`httk.codes.lammps`."""

from httk.core.register import register_code

register_code("lammps", bridge="httk.codes.lammps._bridge", bash_api="httk.codes.lammps:httk-lammps.sh")
