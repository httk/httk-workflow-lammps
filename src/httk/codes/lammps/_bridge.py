"""The ``lammps-*`` subcommands of the private native Bash command bridge.

``httk.workflow._shell_bridge`` mounts these beside its own subcommands through
the ``codes`` registry tier, so each function of ``httk-lammps.sh`` is one
invocation of one command here. A legitimately absent answer returns the
bridge's uniform exit code ``1``; a refused call raises, which the bridge
reports as ``2``. ``lammps-run`` has its own outcome codes, the same as
``vasp-run``: ``0`` completed, ``20`` crashed, ``21`` nonconverged (a
minimization that did not meet its tolerances), ``22`` process failure,
``124`` timeout; ``lammps-diagnose`` exits ``20`` when it found anything, like
``vasp-diagnose``.
"""

import argparse
import json
from pathlib import Path

from httk.workflow.codes import BRIDGE_ABSENT, read_json

from .diagnostics import diagnose_lammps
from .inputs import write_lammps_data
from .outputs import parse_lammps_log
from .reports import run_lammps

_RUN_EXIT = {"completed": 0, "crashed": 20, "nonconverged": 21, "process_failure": 22, "timeout": 124}


def add_commands(commands: "argparse._SubParsersAction[argparse.ArgumentParser]") -> None:
    """Register the ``lammps-*`` subcommands on the bridge's subparsers.

    :param commands: the bridge's subcommand collection.
    """

    run = commands.add_parser("lammps-run")
    run.add_argument("--directory", default=".")
    run.add_argument("--input", default="in.lammps")
    run.add_argument("--log", default="log.lammps")
    run.add_argument("--output", default="lammps.out")
    run.add_argument("--timeout", type=float)
    run.add_argument("argv", nargs=argparse.REMAINDER)
    thermo = commands.add_parser("lammps-thermo")
    thermo.add_argument("--log", default="log.lammps")
    thermo.add_argument("--column", required=True)
    thermo.add_argument("--table", type=int, default=-1)
    diagnose = commands.add_parser("lammps-diagnose")
    diagnose.add_argument("--log", default="log.lammps")
    diagnose.add_argument("--output", default="lammps.out")
    diagnose.add_argument("--json", action="store_true")
    write = commands.add_parser("lammps-write-data")
    write.add_argument("--options", required=True)
    write.add_argument("--data", default="structure.data")


def run_command(arguments: argparse.Namespace) -> int:
    """Run one parsed ``lammps-*`` subcommand.

    :param arguments: the parsed bridge command line.
    :return: the subcommand's exit code.
    """

    command = arguments.command
    if command == "lammps-run":
        argv = arguments.argv[1:] if arguments.argv[:1] == ["--"] else arguments.argv
        if not argv:
            raise ValueError("lammps-run needs the lmp command after --")
        report = run_lammps(
            argv,
            directory=arguments.directory,
            input_file=arguments.input,
            log_file=arguments.log,
            output_file=arguments.output,
            timeout=arguments.timeout,
        )
        print(Path(arguments.directory, "lammps-run-report.json"))
        return _RUN_EXIT[report.classification]
    if command == "lammps-diagnose":
        log = Path(arguments.log)
        diagnostics = diagnose_lammps(log.parent, log=log.name, output=arguments.output)
        if arguments.json:
            print(json.dumps([item.as_mapping() for item in diagnostics], sort_keys=True))
        else:
            for item in diagnostics:
                print(f"{item.code}\t{item.severity}\t{item.summary}")
        return 20 if diagnostics else 0
    if command == "lammps-write-data":
        write_lammps_data(arguments.data, **dict(read_json(Path(arguments.options))))
        return 0
    if command == "lammps-thermo":
        tables = parse_lammps_log(arguments.log).thermo
        try:
            value = tables[arguments.table].last.get(arguments.column)
        except IndexError:
            return BRIDGE_ABSENT
        if value is None:
            return BRIDGE_ABSENT
        # 15 significant digits reproduce every value LAMMPS prints (8 by default) without float noise.
        print(f"{value:.15g}")
        return 0
    raise AssertionError(command)
