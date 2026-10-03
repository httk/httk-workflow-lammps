# *httk-workflow-lammps*

This site documents the *httk-workflow-lammps* module. For the full
documentation of *httk₂*, see [docs.httk.org](https://docs.httk.org).

The module adds LAMMPS support to *httk-workflow*: the Python helpers in
`httk.codes.lammps` (data-file writing, log parsing, diagnostics and supervised
execution), the Bash API a Bash runner sources as
`$HTTK_WORKFLOW_LAMMPS_BASH_API`, and the `lammps-*` bridge commands behind that
API. Installing it registers the `lammps` code with *httk₂* through the
`httk.registry.codes.lammps` registration package. The repository also carries
the example workflow package `lammps.run`.

```{admonition} Quick links
:class: tip

- {doc}`usage` — the Python and Bash API, the example workflow, and the diagnostics
- {doc}`analysis` — joining a LAMMPS dump with its log into canonical-unit samples
- {doc}`reference/index` — the generated API reference
```

## Install

```console
python -m pip install "httk-workflow-lammps[atomistic]"
```

```{toctree}
:maxdepth: 2
:caption: Documentation

usage
analysis
reference/index
```
