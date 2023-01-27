#!/usr/bin/env python3
# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
import sys
import tempfile
from math import ceil
from subprocess import run

from .export import Export
from .logger import logger
from .utils import RUNASTER_ROOT

USAGE = """
    run_sbatch [sbatch-options] FILE.export

This script simply wraps the execution of a study with sbatch:

    sbatch <options from export> .../bin/run_aster FILE.export

'sbatch-options' are passed to 'sbatch' before those deduced from the .export file.
Use 'sbatch --help' for details and example below.
"""

EPILOG = """Example:
    sbatch --wckey=p11yb:aster --partition=bm FILE.export
"""

TEMPLATE = """#!/bin/bash
#SBATCH --job-name={name}

# nbnodes
#SBATCH -N {mpi_nbnodes}

# nb proc MPI total
#SBATCH -n {mpi_nbcpu}

# max walltime
#SBATCH --time="00:00:{time_limit}"

# memory in MB
#SBATCH --mem={memory_node}M

# add `--exclusive` if several nodes
#SBATCH {opt_exclusive}

# redirect output in the current directory
#SBATCH --output=run_aster-%j.txt

{RUNASTER_ROOT}/bin/run_aster {study}
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
    if nbnodes > 1 or cpu_per_node >= 6:
        params["opt_exclusive"] = "--exclusive"


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
    params = dict(
        name=osp.splitext(osp.basename(args.file))[0],
        mpi_nbcpu=export.get("mpi_nbcpu", 1),
        mpi_nbnodes=nbnodes,
        time_limit=export.get("time_limit", 3600),
        memory_limit=export.get("memory_limit", 16384),
        memory_node=export.get("memory_limit", 16384),
        opt_exclusive="",
        study=args.file,
        RUNASTER_ROOT=RUNASTER_ROOT,
    )
    check_parameters(params)

    logger.debug("Parameters: %s", params)
    content = TEMPLATE.format(**params)
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as fobj:
        fobj.write(content)
        script = fobj.name

    logger.info("+ submitted script:\n%s", content)
    try:
        proc = _run(["sbatch"] + sbatch_args + [script])
    finally:
        os.remove(script)
    return proc.returncode


if __name__ == "__main__":
    sys.exit(main())
