# Joining a LAMMPS dump with its log

{py:func}`~httk.codes.lammps.lammps_samples` joins a native atom/custom
text dump with a log parsed by {py:func}`~httk.codes.lammps.parse_lammps_log`.
It uses the exact, streaming `LammpsTrajectory` backend in *httk-atomistic*
(the `atomistic` extra), so the join raises `ImportError` without it. The
yielded samples carry exact structures and canonical-unit thermodynamic values,
which *httk-analyse* consumes without knowing that LAMMPS produced them.

```python
from httk.codes.lammps import lammps_samples

# Use files from the same run and select both sources explicitly.
# samples = lammps_samples(
#     "atoms.lammpstrj", "log.lammps",
#     species={1: "Ar"}, units="metal", dimension=3, segment=0, table_index=0,
#     columns=["Temp", "PotEng", "TotEng", "Volume"],
#     observables=["atom_ids", "unwrapped_positions", "velocities", "time"],
#     trajectory_options={"timestep": "0.001"}, join="strict",
# )
# for sample in samples:
#     print(sample.step, sample.thermo)
```

Print
`id type xu yu zu vx vy vz` in the dump and `step atoms temp pe etotal vol`
in the thermo table. `dump_modify ... time yes` can supply actual times;
otherwise provide a constant timestep only when it represents that run.
Changing timesteps require printed times. Integer atom types need an explicit
species mapping. Source decimal data remain exact in structures and trajectory
observables; thermo conversion uses the log parser's float values.

## Matching and conversion

`intersection` emits shared steps in dump order. `strict` requires identical
ordered step sequences and raises for a cadence mismatch. Neither interpolates
or matches rows by position. A strict iterator may yield earlier samples before
it discovers a later mismatch; consume the iterator fully before accepting a
complete run. A dump segment begins on a non-increasing step. Select the matching
log table explicitly, including across restarted runs. The user is responsible
for selecting the same physical run; step equality alone cannot prove this.

Supply `dimension=3` explicitly for these bulk quantities. Two-dimensional
LAMMPS pressure and volume use area conventions and are rejected by this join.

Canonical outputs use eV, angstrom, ps and K. Positive pressure is compression;
converted stress is tensile-positive. For stress components, explicitly set
`stress_basis_matches_dump=True` only after verifying the Cartesian bases:
general triclinic output can rotate the dump relative to default thermo tensors.
This join does not infer that rotation.

For normalized energies, the global `Atoms` thermo column or an explicit
`total_atom_count` is required. A dump may include only a group, so its atom
count cannot establish the normalization. Volume is printed as the full box
volume even when energies are normalized. LJ units need positive explicit
length/time/energy scales in `trajectory_options['lj_scales']`. Custom variables
need an explicit `ThermoConversion(name, unit, scale, per_atom)`; their meaning
and normalization are never inferred from column names. The `per_atom` flag
means that the printed value must be multiplied by the global atom count.

See the official [dump format](https://docs.lammps.org/dump.html),
[thermo output](https://docs.lammps.org/thermo_style.html), and
[unit conventions](https://docs.lammps.org/units.html).
