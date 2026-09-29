# Test data

Captured LAMMPS (22 Jul 2025 - Update 6) runs, serial with
`OMP_NUM_THREADS=1`, run as `lmp -in <input> -log <log>`. The `.log` files are
the `-log` output, the `.out` files the captured standard output.

| File | What it is |
| --- | --- |
| `lj.in`, `lj.log`, `lj.out` | Lennard-Jones fcc melt, 108 atoms, 100 NVE steps; final thermo row `100 0.57389677 -6.1447617 -5.2918874 -2.0049305` (`Step Temp PotEng TotEng Press`), `Total wall time:`, exit status 0 |
| `two_runs.in`, `two_runs.log` | the melt as two `run` commands with different `thermo_style`s; `fix print` lines interleave with the thermo rows, and the second `thermo_style` prints a `WARNING:` |
| `min.in`, `min.log` | `minimize` of the perfect lattice: `Stopping criterion = energy tolerance` |
| `min_maxiter.in`, `min_maxiter.log` | `minimize` of a randomly displaced lattice capped at 5 iterations: `Stopping criterion = max iterations` |
| `err.in`, `err.log`, `err.out` | an unknown command: `ERROR: Unknown command: this_is_not_a_lammps_command foo bar (src/input.cpp:315)`, exit status 1 |
| `noinput.out` | standard output of `lmp -in nothere`: `ERROR on proc 0: Cannot open input script nothere: ...`, an empty log, exit status 1 |

The inputs are written for these tests; the logs and outputs are program output
of those runs.
