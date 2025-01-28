# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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
:py:mod:`toolbox` --- Toolbox for the developers
------------------------------------------------
"""

import pathlib
from .config import CFG
from .utils import run_command
import os

def make_shared(lib, src, *args):
    """Build a shared library from a fortran source file.

    Arguments:
        lib (str): Result library.
        src (str): Fortran source file.
        args (list[str]): Optional arguments passed to the linker.

    Returns:
        int: exit code.
    """
    fc = CFG.get("FC")
    fcp = pathlib.Path(fc)
    if not fcp.exists():
        fc_env = os.environ.get("FC")
        if fc_env is None:
            raise FileNotFoundError(f"Fortran compiler not found: {fc}")
        fcp = pathlib.Path(fc_env)
    cmd = [fcp.as_posix()]
    cmd.extend(CFG.get("FCFLAGS"))
    cmd.extend(["-shared", "-o", lib, src])
    cmd.extend(args)
    print("INFO make_shared command line:", " ".join(cmd))
    return run_command(cmd)
