#!/usr/bin/env python3
# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

"""
``bin/run_aster`` --- Script to execute code_aster from a ``.export`` file
--------------------------------------------------------------------------

``bin/run_aster`` executes a code_aster study from the command line.
The parameters and files used by the study are defined in a ``.export`` file
(see :py:mod:`~run_aster.export` for description of the syntax of the
``.export``).

For parallel executions, these two forms are equivalent:

.. code-block:: sh

    bin/run_aster path/to/file.export

or:

.. code-block:: sh

    mpiexec -n 4 bin/run_aster path/to/file.export

Using the first syntax, ``bin/run_aster`` re-runs with ``mpiexec`` itself using
the second syntax (``mpiexec`` syntax is provided by the configuration, see
:py:mod:`~run_aster.config`).

``bin/run_aster`` only runs its own version, those installed at the level of
the ``bin`` directory; unlike ``as_run`` where the same instance of ``as_run``
executes several versions of code_aster.
This makes ``bin/run_aster`` simpler and allows per version settings
(see :py:mod:`~run_aster.config` for more informations about the configuration
of a version).

More options are available to execute code_aster with an interactive Python
interpreter, to prepare a working directory and to start manually (through
a debugger for example)...
For example, executing ``bin/run_aster`` with no ``.export`` file starts an
interactive Python interpreter.

See ``bin/run_aster --help`` for the available options.

"""

# aslint: disable=C4009
# C4009: in a string, imported here

import argparse
import os
import os.path as osp
import shutil
import sys
import tempfile
from subprocess import run

from .command_files import AUTO_IMPORT
from .config import CFG
from .export import Export, File, split_export
from .logger import DEBUG, WARNING, logger
from .run import RunAster, create_temporary_dir, get_procid
from .status import StateOptions, Status
from .utils import RUNASTER_ROOT

try:
    import debugpy

    HAS_DEBUGPY = True
except ImportError:
    HAS_DEBUGPY = False

USAGE = """
    run_aster [options] [EXPORT]

"""

EPILOG = """
The time limit is automatically increased to 24 hours for "interactive" usages
as '--interactive', '--env', '--no-comm'.
"""

# for interactive executions (using IPython)
SAVE_ARGV = """
import sys

from code_aster.Utilities import ExecutionParameter
ExecutionParameter().set_argv(sys.argv)
"""


def parse_args(argv):
    """Parse command line arguments.

    Arguments:
        argv (list): List of command line arguments.
    """
    # command arguments parser
    parser = argparse.ArgumentParser(
        usage=USAGE,
        epilog=EPILOG,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version", action="store_true", help="show code_aster version"
    )
    parser.add_argument(
        "-g",
        "--debug",
        action="store_true",
        help="print debugging information (same as " "DEBUG=1 in environment)",
    )
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        help="inspect interactively after running script "
        "instead of calling FIN command",
    )
    parser.add_argument(
        "--env",
        action="store_true",
        help="do not execute, only prepare the working "
        "directory ('--wrkdir' is required)",
    )
    parser.add_argument(
        "-w", "--wrkdir", action="store", help="use this directory as working directory"
    )
    parser.add_argument(
        "--all-procs",
        dest="only_proc0",
        action="store_false",
        default=None,
        help="show all processors output",
    )
    parser.add_argument(
        "--only-proc0",
        dest="only_proc0",
        action="store_true",
        default=None,
        help="only processor #0 is writing on stdout",
    )
    parser.add_argument(
        "--no-mpi",
        dest="auto_mpiexec",
        action="store_false",
        help="if '%(prog)s' is executed with a parallel "
        "version but not under 'mpiexec', it is "
        "automatically restart with "
        "'mpiexec -n N %(prog)s ...'; use '--no-mpi' "
        "to not do it",
    )
    parser.add_argument(
        "-t", "--test", action="store_true", help="execution of a testcase"
    )
    parser.add_argument(
        "--ctest",
        action="store_true",
        help="testcase execution inside ctest (implies "
        "'--test'), the 'mess' file is saved into "
        "the current directory (which is '--resutest' "
        "directory for 'run_ctest') and is not duplicated "
        "on stdout.",
    )
    parser.add_argument(
        "--time_limit",
        dest="time_limit",
        type=float,
        action="store",
        default=None,
        help="override the time limit (may also be changed by "
        "FACMTPS environment variable)",
    )
    parser.add_argument(
        "--memory_limit",
        dest="memory_limit",
        type=float,
        action="store",
        help="override the memory limit in MB",
    )
    parser.add_argument(
        "--no-comm",
        action="store_true",
        help="do not execute the `.comm` files but start an "
        "interactive Python session. Execute "
        "`code_aster.init()` to copy data files.",
    )
    parser.add_argument(
        "--exectool",
        action="store",
        help="wrap code_aster execution using this tool "
        "(debugger, valgrind, custom command...)",
    )
    parser.add_argument(
        "--debugpy-runner", action="store", type=int, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--debugpy-rank", action="store", type=int, default=0, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--status-file", action="store", dest="statusfile", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "export",
        metavar="EXPORT",
        nargs="?",
        help="Export file defining the calculation. "
        "Without file, it starts an interactive Python "
        "session.",
    )

    args = parser.parse_args(argv)
    if args.ctest:
        logger.setLevel(WARNING)
    if args.debug:
        logger.setLevel(DEBUG)
        os.environ["DEBUG"] = str(os.environ.get("DEBUG") or 1)
    if args.version:
        tag = CFG.get("version_tag")
        sha1 = CFG.get("version_sha1")[:12]
        logger.info(f"code_aster {tag} ({sha1})")
        parser.exit(0)
    if args.env and not args.wrkdir:
        parser.error("Argument '--wrkdir' is required if '--env' is enabled")
    if args.debugpy_runner and not HAS_DEBUGPY:
        parser.error("can not import 'debugpy'")
    return args


def main(argv=None):
    """Entry point for code_aster runner.

    Arguments:
        argv (list): List of command line arguments.
    """
    argv = argv or sys.argv[1:]
    args = parse_args(argv)

    procid = 0
    if CFG.get("parallel", 0):
        procid = get_procid()

    export = Export(args.export, " ", test=args.test or args.ctest, check=False)
    make_env = args.env or "make_env" in export.get("actions", [])
    need_split = len(export.commfiles) > 1
    if need_split and (CFG.get("parallel", 0) and procid >= 0):
        logger.error(
            "Can not execute several comm files under MPI runner. "
            "Let run_aster split the export file or change the export file."
        )
    need_mpiexec = procid < 0 and args.auto_mpiexec
    logger.debug(
        "parallel: {0}, procid: {1}".format(CFG.get("parallel", False), procid)
    )
    logger.debug("nbcomm: {0}".format(len(export.commfiles)))
    if (
        args.debugpy_runner
        and procid == args.debugpy_rank
        and not (need_split or need_mpiexec)
    ):
        debugpy.listen(("localhost", args.debugpy_runner))
        print("Waiting for debugger attach")
        debugpy.wait_for_client()
        debugpy.breakpoint()

    tmpf = None
    if args.no_comm:
        for comm in export.commfiles:
            export.remove_file(comm)
    if not args.export or args.no_comm:
        args.interactive = True
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as fobj:
            fobj.write(AUTO_IMPORT.format(starter=SAVE_ARGV))
            export.add_file(File(fobj.name, filetype="comm", unit=1))
            tmpf = fobj.name
    if args.ctest:
        args.test = True
        basename = osp.splitext(osp.basename(args.export))[0]
        add = {6: "mess", 15:"code"}
        for unit, typ in add.items():
            export.add_file(File(osp.abspath(basename + "." + typ),
                                 filetype=typ,
                                 unit=unit,
                                 resu=True))
    if args.time_limit:
        export.set_time_limit(args.time_limit)
    # use FACMTPS from environment
    try:
        mult = float(os.environ.get("FACMTPS", 1))
        limit = export.get("time_limit", 86400.0) * mult
        export.set_time_limit(limit)
    except ValueError:
        pass
    if args.interactive:
        export.set("interact", True)
    if args.interactive or make_env:
        export.set_time_limit(86400.0)
    if args.memory_limit:
        export.set_memory_limit(args.memory_limit)
    if export.get("no-mpi"):
        args.auto_mpiexec = False
    export.check()

    if args.only_proc0 is None:
        args.only_proc0 = CFG.get("only-proc0", False)

    wrkdir = args.wrkdir or create_temporary_dir(dir=CFG.get("tmpdir"))
    try:
        if need_split or need_mpiexec:
            run_aster = osp.join(RUNASTER_ROOT, "bin", "run_aster")
            expdir = create_temporary_dir(
                dir=os.getenv("HOME", "/tmp") + "/.tmp_run_aster"
            )
            statfile = osp.join(expdir, "__status__")
            for exp_i in split_export(export):
                fexp = osp.join(expdir, "export." + str(exp_i.get("step")))
                exp_i.write_to(fexp)
                argv_i = [i for i in argv if i != args.export]
                if not args.wrkdir:
                    argv_i.append("--wrkdir")
                    argv_i.append(wrkdir)
                argv_i.extend(["--status-file", statfile])
                argv_i.append(fexp)
                cmd = f"{run_aster} {' '.join(argv_i)}"
                if need_mpiexec:
                    args_cmd = dict(mpi_nbcpu=export.get("mpi_nbcpu", 1), program=cmd)
                    cmd = CFG.get("mpiexec").format(**args_cmd)
                logger.info("Running: " + cmd)
                proc = run(cmd, shell=True)
                status = Status.load(statfile)
                if proc.returncode != 0 and not status.is_completed():
                    break
            shutil.rmtree(expdir)
            return proc.returncode

        if args.only_proc0 and procid > 0:
            logger.setLevel(WARNING)
        opts = {}
        opts["test"] = args.test
        opts["env"] = make_env
        opts["tee"] = not args.ctest and (not args.only_proc0 or procid == 0)
        opts["interactive"] = args.interactive
        if args.exectool:
            wrapper = CFG.get("exectool", {}).get(args.exectool)
            if not wrapper:
                logger.warning(
                    f"'{args.exectool}' is not defined in your "
                    f"configuration, it is used as a command line."
                )
                wrapper = args.exectool
            opts["exectool"] = wrapper
        calc = RunAster.factory(export, **opts)
        status = calc.execute(wrkdir)
        if args.statusfile:
            status.save(args.statusfile)
        if tmpf and not opts["env"]:
            os.remove(tmpf)
    finally:
        if not args.wrkdir:
            shutil.rmtree(wrkdir)
    return status.exitcode


if __name__ == "__main__":
    sys.exit(main())
