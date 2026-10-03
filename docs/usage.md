# Using the LAMMPS helpers

*httk-workflow-lammps* ships the LAMMPS helpers that workflow runners are built
on, in two languages: the Python package {py:mod}`httk.codes.lammps` and the
Bash LAMMPS API, whose `httk_lammps_*` functions call the same code through the
*httk-workflow* shell bridge.

## Install

```console
python -m pip install "httk-workflow-lammps[atomistic]"
```

The distribution depends on *httk-core* and *httk-workflow*. Installing it
registers the `lammps` code through the `httk.registry.codes.lammps`
registration package, which makes the `lammps-*` bridge commands and the Bash
API available to every job the manager starts; nothing needs to be configured.
Writing data files from structures needs *httk-atomistic*, the `atomistic`
extra; parsing, diagnostics and running do not.

## Python

```python
from httk.codes.lammps import run_lammps, write_lammps_data

write_lammps_data("structure.data", structure="POSCAR", masses={"Ar": 39.948})
report = run_lammps(["mpirun", "-np", "4", "lmp"], timeout=3600)
if report.ok:
    print(report.result.thermo[-1].last["TotEng"])
else:
    print(report.classification, [item.code for item in report.diagnostics])
```

- {py:func}`~httk.codes.lammps.write_lammps_data` writes a LAMMPS data file for
  `read_data` from a POSCAR/CIF path or an *httk* structure: the cell rotated
  into the LAMMPS restricted-triclinic frame (an `xy xz yz` tilt line unless
  the cell is orthogonal), atom types numbered by first appearance, Cartesian
  positions in Å, and a `Masses` section when every species has a mass (from
  `masses` or the structure). The values are floats of the structure's exact
  ones. Only `atom_style atomic` is written.
- {py:func}`~httk.codes.lammps.parse_lammps_log` returns a
  {py:class}`~httk.codes.lammps.LammpsResult`: one
  {py:class}`~httk.codes.lammps.ThermoTable` per `run` or `minimize` command
  (its header columns, float rows, and `last` row as a mapping; lines that
  interleave with the rows are skipped), whether the `Total wall time:` end was
  reached, the `ERROR` and `WARNING` lines, and the stopping criterion of the
  last minimization. Only the tabular thermo styles (`one`, `custom`) are read
  as tables.
- {py:func}`~httk.codes.lammps.run_lammps` runs the command with
  `-in in.lammps -log log.lammps` under the *httk-workflow* process supervisor,
  saves standard output as `lammps.out`, and returns a
  {py:class}`~httk.codes.lammps.LammpsRunReport` classified as `completed`,
  `crashed`, `nonconverged`, `process_failure` or `timeout`, also written to
  `lammps-run-report.json`.
- {py:func}`~httk.codes.lammps.diagnose_lammps` diagnoses a finished run from
  its log, plus the error lines of its standard output (an unreadable input
  script is reported only there).

## Diagnostics

| Code | Severity | Meaning |
| --- | --- | --- |
| `lammps.error` | fatal | LAMMPS stopped with an `ERROR:` or `ERROR on proc N:` line; the summary is the line |
| `lammps.incomplete` | error | no `Total wall time:` line and no error, e.g. a killed process |
| `lammps.minimization_not_converged` | error | the last minimization stopped on `max iterations` or `max force evaluations` |

`WARNING:` lines are not diagnostics; they are kept in the result's `warnings`.
A clean run has no diagnostics.

## Bash

The manager exports the path of the LAMMPS API as
`HTTK_WORKFLOW_LAMMPS_BASH_API` when *httk-workflow-lammps* is installed, so a
Bash runner guards it and sources it after the generic library:

```bash
source "$HTTK_WORKFLOW_BASH_API"
: "${HTTK_WORKFLOW_LAMMPS_BASH_API:?install httk-workflow-lammps}"
source "$HTTK_WORKFLOW_LAMMPS_BASH_API"

httk_lammps_write_data --options options.json   # the write_lammps_data keywords as JSON
httk_lammps_run --timeout 3600 -- mpirun -np 4 lmp
energy=$(httk_lammps_thermo --column TotEng)
```

| Function | Bridge command | Exit status |
| --- | --- | --- |
| `httk_lammps_write_data --options FILE [--data structure.data]` | `lammps-write-data` | `0` |
| `httk_lammps_run [--directory] [--input] [--log] [--output] [--timeout] -- CMD...` | `lammps-run` | `0` completed, `20` crashed, `21` nonconverged (minimization), `22` process failure, `124` timeout (as `vasp-run`); prints the report path |
| `httk_lammps_thermo --column NAME [--table N] [--log log.lammps]` | `lammps-thermo` | `0` and the last value of the column in table `N` (default `-1`, the last), `1` when there is none |
| `httk_lammps_diagnose [--log log.lammps] [--output lammps.out] [--json]` | `lammps-diagnose` | `0` clean, `20` when it printed diagnostics |

A refused call (for example a missing log file) exits `2`. The API sets
`HTTK_LAMMPS_BASH_API_VERSION=1`.

## The example workflow

The repository's `workflows/lammps-run` is the workflow package `lammps.run`:
one Python runner step that stages the `script` input as `in.lammps`, writes the
optional `structure` input as `structure.data` (for a script that says
`read_data structure.data`), runs LAMMPS, and fails with the first diagnostic
code when the run is not clean. Install it with `httk plugin install` of the
repository, or use it directly with `--workflow-dir`:

```console
httk workspace settings set --key lammps.command --value 'mpirun -np 4 lmp' WORKSPACE
httk job new --workflow lammps.run --input script=in.lammps --input structure=POSCAR
httk workflow run
httk collect --into results.sqlite
```

The setting `lammps.command` (default `lmp`) says how to run LAMMPS. The data
file carries no masses unless the structure does, so the script sets them with
`mass`.

## Collecting

`lammps.run` declares no outputs and has no collector, so collecting it stores
only the job's run: a run-only provenance record naming the workflow and its
inputs, with the files left in the job's workdir.

### Recognized calculations

A finished LAMMPS run that was not started by a workspace is collected by the
registered `lammps.calculation` collector:
`httk.workflow.collect_tree(root)` finds every directory holding exactly one
log (any name, compressed or not) whose first five lines carry the
`LAMMPS (` banner, together with exactly one input script, and collects the
output `average_total_energy`. The identity is a digest of the input script and the files it reads,
so moving the directory keeps it.

- **Identity.** The digest covers the script and the plain relative files in
  the directory named by its `read_data`, `read_restart` and `include` lines
  (first argument), and by any argument of its `pair_coeff` (potential files)
  and `molecule` lines (names containing `$` variables or `/` are ignored). Runs that
  differ only by command-line `-var` values are indistinguishable by design.
- **Screen copies.** LAMMPS prints the banner on screen too: when a captured
  `slurm-N.out` sits beside the log, the file that echoes the input (else
  `log.*`/`*.log`) is taken as the log.

- **Averaging rule.** LAMMPS prints no average, so the value is the arithmetic
  mean, computed by *httk₂*, of the printed `TotEng` column over the rows of
  the **last `run`** thermo table; `minimize` tables and a trailing `run 0`
  are skipped. A staged `run 500` then `run 500` therefore averages only the
  last segment. The units and `thermo_modify norm` state used are those in
  force at that table (a later `units` or `clear` is followed); the claim
  decision reads them from the first 2000 echoed lines only, and the collect
  step degrades if the full log disagrees. The
  `thermo_style` must print `etotal` (`thermo_style multi` and `yaml` are not
  parsed). A log without a completed run (`Total wall time:`), without a `run`
  table or without `TotEng` is claimed and degraded.
- **Units and normalization.** The `units` style is read from the echoed input
  at the top of the log. LAMMPS defaults thermo normalization to yes for `lj`
  and no for its other known styles; each `thermo_style` resets it to that
  default, and an explicit `thermo_modify norm` overrides it until the next
  reset. `clear`, missing setup, an unrecognized style, or an invalid echoed
  value leaves the state unknown. The `average_total_energy_ev` helper refuses
  unknown or per-atom state because the printed `TotEng` may not be a system
  total. `metal` (eV) and `real` (kcal/mol, converted with
  {py:data}`~httk.codes.lammps.KCAL_MOL_TO_EV`, the thermochemical kcal times
  the CODATA kJ/mol to eV factor) are accepted; every other style, or a log
  with `echo none`, is reported as unclaimed, as is an input with
  `thermo_modify norm yes` (per-atom energies, not the whole system).
- **File names.** The collector is selected by a `*.lammps` file in the
  directory, so the workflow's `in.lammps`/`log.lammps` pair is found; the log
  is then located by its banner, whatever its name. The input script is the
  single `in.*` or `*.in` file beside the log. A directory whose only
  `.lammps` file is absent (for example `md.in` with `run.log`) is not
  recognized, and several logs or several scripts are reported as unclaimed.
  {py:func}`~httk.codes.lammps.collect.find_outputs` is the banner-based finder.
