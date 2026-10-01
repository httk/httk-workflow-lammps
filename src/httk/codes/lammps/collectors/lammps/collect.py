"""Collect hook for the ``lammps.calculation`` collector.

The directory holds one LAMMPS log, found by its banner, whatever its name.
"""

from httk.workflow.collecting import JobRecord

from httk.codes.lammps.collect import find_outputs, read_average_total_energy


def collect(record: JobRecord):
    """Return the mean total energy of the calculation.

    :param record: The stand-in job record of the directory.
    :return: The ``average_total_energy`` output role.
    """
    if record.workdir is None:
        raise ValueError("the job has no working directory")
    (output,) = find_outputs(record.workdir)
    return {"average_total_energy": read_average_total_energy(output)}
