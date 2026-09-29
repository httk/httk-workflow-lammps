#!/usr/bin/env bash

# Native httk LAMMPS Bash API, version 1. Source httk-workflow.sh first.
#
# Every function is one lammps-* bridge subcommand, and every option of that
# subcommand is available here: the arguments are passed through untouched.
#
#   httk_lammps_write_data --options OPTIONS.json [--data structure.data]
#   httk_lammps_run [--directory .] [--input in.lammps] [--log log.lammps] [--output lammps.out]
#                   [--timeout S] -- lmp ...
#       prints the report path; exits 0 completed, 20 crashed, 21 nonconverged
#       (minimization), 22 process failure, 124 timeout
#   httk_lammps_thermo --column NAME [--table N] [--log log.lammps]
#       prints the last value of the column in thermo table N (default -1, the
#       last; negative counts from the end); exits 1 when there is none
#   httk_lammps_diagnose [--log log.lammps] [--output lammps.out] [--json]
#       exits 20 when it found anything
HTTK_LAMMPS_BASH_API_VERSION=1

_httk_lammps_require_workflow_api() {
    if ! declare -F _httk_workflow_bridge >/dev/null 2>&1; then
        printf 'httk-workflow: source HTTK_WORKFLOW_BASH_API before HTTK_WORKFLOW_LAMMPS_BASH_API\n' >&2
        return 2
    fi
}

httk_lammps_write_data() {
    _httk_lammps_require_workflow_api || return
    _httk_workflow_bridge lammps-write-data "$@"
}

httk_lammps_run() {
    _httk_lammps_require_workflow_api || return
    _httk_workflow_bridge lammps-run "$@"
}

httk_lammps_thermo() {
    _httk_lammps_require_workflow_api || return
    _httk_workflow_bridge lammps-thermo "$@"
}

httk_lammps_diagnose() {
    _httk_lammps_require_workflow_api || return
    _httk_workflow_bridge lammps-diagnose "$@"
}
