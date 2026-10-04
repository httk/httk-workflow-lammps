# Joining a LAMMPS dump with its log

{py:func}`~httk.codes.lammps.lammps_samples` joins a native atom/custom
text dump with a log parsed by {py:func}`~httk.codes.lammps.parse_lammps_log`.
It uses the exact, streaming `LammpsTrajectory` backend in *httk-atomistic*
(the `atomistic` extra), so the join raises `ImportError` without it. The
yielded samples carry exact structures and thermodynamic values in the units fixed by their OPTIMADE property definitions,
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

Each `sample.thermo` entry is a pair of a property definition name
(`total_energy`, `potential_energy`, `kinetic_energy`, `enthalpy`, `temperature`,
`volume`, `pressure`, `stress_tensor`) and a bare value in the unit of that
definition: eV for whole-cell energies, K, angstrom^3, and GPa. No unit strings
are carried; the conversion happens once, here. Pairs follow the order of
`columns`. Positive `pressure` is compression. The six stress columns
`Pxx, Pyy, Pzz, Pyz, Pxz, Pxy` must be selected together and yield a single
`stress_tensor` value, a tensile-positive Voigt tuple `(xx, yy, zz, yz, xz, xy)`
in GPa, placed at the position of the first of them. A custom conversion cannot produce a stress tensor. For stress, explicitly set
`stress_basis_matches_dump=True` only after verifying the Cartesian bases:
general triclinic output can rotate the dump relative to default thermo tensors.
This join does not infer that rotation.

For normalized energies, the global `Atoms` thermo column or an explicit
`total_atom_count` is required. A dump may include only a group, so its atom
count cannot establish the normalization. Volume is printed as the full box
volume even when energies are normalized. LJ units need positive explicit
length/time/energy scales in `trajectory_options['lj_scales']`. Custom variables
need an explicit `ThermoConversion(definition=<property definition IRI>,
scale=<factor>, per_atom=<bool>)`; the caller chooses `scale` so that the printed
value times it is in the definition's unit. Meaning and normalization are never
inferred from column names. The `per_atom` flag
means that the printed value must be multiplied by the global atom count.

See the official [dump format](https://docs.lammps.org/dump.html),
[thermo output](https://docs.lammps.org/thermo_style.html), and
[unit conventions](https://docs.lammps.org/units.html).
