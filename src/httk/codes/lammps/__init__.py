"""LAMMPS (``lmp``) support for *httk₂* workflows: the *httk-workflow-lammps* package.

``inputs`` writes LAMMPS data files, ``outputs`` parses the LAMMPS log,
``diagnostics`` classifies a finished run, and ``reports`` runs it under
supervision. This package is a thin facade re-exporting their surface. The
example workflow package ``workflows/lammps-run`` in this distribution's
repository builds on it.
"""

from httk.core import register_citation

register_citation(
    applies_to="Calculations with LAMMPS",
    references=(
        {
            "authors": ({"name": "Aidan P. Thompson"},),
            "note": "Thompson et al.; the DOI record lists every author",
            "title": (
                "LAMMPS - a flexible simulation tool for particle-based materials modeling"
                " at the atomic, meso, and continuum scales"
            ),
            "journal": "Computer Physics Communications",
            "volume": "271",
            "pages": "108171",
            "year": "2022",
            "doi": "10.1016/j.cpc.2021.108171",
            "bib_type": "article",
        },
    ),
)

from .diagnostics import diagnose_lammps
from .inputs import write_lammps_data
from .outputs import LammpsResult, ThermoTable, parse_lammps_log
from .reports import LammpsRunReport, run_lammps

__all__ = [
    "LammpsResult",
    "LammpsRunReport",
    "ThermoTable",
    "diagnose_lammps",
    "parse_lammps_log",
    "run_lammps",
    "write_lammps_data",
]
