# httk-workflow-lammps

![Status: Early beta](https://img.shields.io/badge/status-early--beta-orange)

> **⚠️ EARLY BETA**
>
> This is an early beta release of *httk₂*. The organization of the packages
> and their APIs should not yet be regarded as stable, and may change between
> releases.

*httk-workflow-lammps* adds LAMMPS support to
[*httk-workflow*](https://github.com/httk/httk-workflow), the workflow engine of
[*httk₂*](https://github.com/httk/httk2). It provides `httk.codes.lammps`:
writing LAMMPS data files from structures, parsing the LAMMPS log (thermo
tables, errors, warnings, minimization stops), stable diagnostics, and
supervised execution with a classified run report; and the Bash API that
exposes the same helpers to Bash runners. Installing it registers the `lammps`
code with *httk₂*; nothing needs to be configured.

## Install

```console
python -m pip install "httk-workflow-lammps[atomistic]"
```

The `atomistic` extra (*httk-atomistic*) is needed only to write data files
from structures.

## Use

In a Python runner:

```python
from httk.codes.lammps import run_lammps, write_lammps_data

write_lammps_data("structure.data", structure="POSCAR", masses={"Ar": 39.948})
report = run_lammps(["lmp"])  # runs lmp -in in.lammps -log log.lammps
print(report.classification, report.result.thermo[-1].last)
```

In a Bash runner, whose manager exports the path of the LAMMPS API:

```bash
source "$HTTK_WORKFLOW_BASH_API"
: "${HTTK_WORKFLOW_LAMMPS_BASH_API:?install httk-workflow-lammps}"
source "$HTTK_WORKFLOW_LAMMPS_BASH_API"
httk_lammps_run -- lmp
temperature=$(httk_lammps_thermo --column Temp)
```

A complete example workflow package, `lammps.run`, is in
[`workflows/lammps-run`](workflows/lammps-run); `httk plugin install` of this
repository installs it. The API is documented in [`docs/usage.md`](docs/usage.md)
and at [docs.httk.org/httk-workflow-lammps](https://docs.httk.org/httk-workflow-lammps/).

## Running tests

`make test` runs the normal profile; `make ci` runs formatting, lint, both type
checkers and the extended tests. The tests that run the real LAMMPS do so only
when `HTTK_TEST_LAMMPS_COMMAND` names it or `lmp` is on `PATH`, and skip
otherwise.
