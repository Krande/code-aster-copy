#!/usr/bin/env python3
# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
``bin/run_sbatch`` --- Script to execute code_aster using ``sbatch``
--------------------------------------------------------------------

``bin/run_sbatch`` executes code_aster using ``sbatch``.

Usage:

.. code-block:: sh

    bin/run_sbatch [sbatch-options] FILE.export

`sbatch-options` are passed to ``sbatch``.

See ``bin/run_sbatch --help`` for the available options.

"""

import argparse
import os
import os.path as osp
import re
import stat
import sys
import tempfile
from math import ceil
from subprocess import run

from .export import Export
from .logger import logger
from .utils import RUNASTER_ROOT
from .config import CFG

USAGE = """
    run_sbatch [sbatch-options] FILE.export

This script simply wraps the execution of a study with sbatch:

    sbatch <options from export> .../bin/run_aster FILE.export

'sbatch-options' are passed to 'sbatch' before those deduced from the .export file.
Use 'sbatch --help' for details and example below.
"""

EPILOG = """Example:
    run_sbatch --wckey=p11yb:aster --partition=bm FILE.export
or:
    export SBATCH_WCKEY=p11yb:aster
    export SBATCH_PARTITION=bm
    run_sbatch FILE.export
"""

TEMPLATE = """#!/bin/bash
#SBATCH --job-name={name}

# number of nodes
#SBATCH --nodes={mpi_nbnodes}

# number of MPI processes
#SBATCH --ntasks={mpi_nbcpu}

# number of threads per MPI process
#SBATCH --cpus-per-task={nbthreads} --threads-per-core=1

# max walltime
#SBATCH --time="00:00:{time_limit}"

# memory in MB
#SBATCH --mem={memory_node}M

# add `--exclusive` if several nodes, define `--partition=...`
#SBATCH {options}

# redirect output in the current directory
#SBATCH --output={output}

{RUNASTER_ROOT}/bin/run_aster {run_aster_options} {study}
"""


def parse_args(argv):
    """Parse command line arguments.

    Arguments:
        argv (list): List of command line arguments.
    """
    parser = argparse.ArgumentParser(
        usage=USAGE, epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-n", "--dry-run", action="store_true", help="do not execute, just show the script content"
    )
    parser.add_argument(
        "--output", action="store", help="output file (default: <export filename>-%%j.txt)"
    )
    parser.add_argument(
        "--run_aster_option",
        dest="opts",
        action="append",
        default=[],
        help="option to be passed to run_aster, can be repeated "
        "(example: --run_aster_option='--only-proc0')",
    )
    parser.add_argument(
        "--ctest",
        dest="opts",
        action="append_const",
        const="--ctest",
        help="shortcut for --run_aster_option='--ctest'",
    )
    parser.add_argument(
        "--time_limit",
        dest="time_limit",
        type=float,
        action="store",
        default=None,
        help="override the time limit in seconds",
    )
    parser.add_argument(
        "--memory_limit",
        dest="memory_limit",
        type=float,
        action="store",
        default=None,
        help="override the memory limit in MB",
    )
    parser.add_argument(
        "file", metavar="FILE.export", help="Export file (.export) defining the calculation."
    )

    args, others = parser.parse_known_args(argv)
    return args, others


def _run(cmd):
    logger.debug("execute: %s", " ".join(cmd))
    return run(cmd)


def check_parameters(params):
    """Check parameters consistency.

    Arguments:
        params (dict): Current parameters.
    """
    nbnodes = params["mpi_nbnodes"]
    cpu_per_node = ceil(params["mpi_nbcpu"] / nbnodes)
    params["memory_node"] = int(cpu_per_node * params["memory_limit"])
    params["time_limit"] = int(params["time_limit"])
    if nbnodes > 1 or cpu_per_node >= 6 or "performance" in params["testlist"]:
        params["options"] += " --exclusive"
    if "bm" in params["testlist"]:
        params["options"] += " --partition=bm"


def main(argv=None):
    """Entry point for sbatch wrapper.

    Arguments:
        argv (list): List of command line arguments.
    """
    args, sbatch_args = parse_args(argv or sys.argv[1:])

    export = Export(args.file)

    # not anymore in Export
    with open(args.file, "r") as fobj:
        text = fobj.read()
    re_nod = re.compile("P +mpi_nbnoeud +([0-9]+)", re.M)
    nbnodes = 1
    mat = re_nod.search(text)
    if mat:
        nbnodes = int(mat.group(1))

    # initialized with default values
    addmem = CFG.get("addmem", 0.0)
    memory = args.memory_limit or export.get("memory_limit", 16384)
    memory += addmem
    if args.time_limit:
        args.opts.append(f"--time_limit={args.time_limit}")
    if args.memory_limit:
        args.opts.append(f"--memory_limit={args.memory_limit}")
    params = {
        "name": osp.splitext(osp.basename(args.file))[0],
        "mpi_nbcpu": export.get("mpi_nbcpu", 1),
        "mpi_nbnodes": nbnodes,
        "nbthreads": export.get("ncpus", 1),
        "time_limit": args.time_limit or export.get("time_limit", 3600),
        "memory_limit": memory,
        "memory_node": None,
        "options": "",
        "study": args.file,
        "run_aster_options": " ".join(args.opts),
        "RUNASTER_ROOT": RUNASTER_ROOT,
        "testlist": export.get("testlist", []),
    }
    params["output"] = args.output or params["name"] + "-%j.txt"
    check_parameters(params)

    logger.debug("Parameters: %s", params)
    content = TEMPLATE.format(**params)
    with tempfile.NamedTemporaryFile(prefix="batch", suffix=".sh", mode="w", delete=False) as fobj:
        fobj.write(content)
        script = fobj.name

    logger.info("+ submitted script:\n%s", content)
    os.chmod(script, stat.S_IRWXU)
    if args.dry_run:
        logger.info("+ filename: %s", script)
        return 0
    try:
        proc = _run(["sbatch"] + sbatch_args + [script])
    finally:
        os.remove(script)
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main())
