# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

import code_aster
from code_aster.Commands import *

from mpi4py import MPI

code_aster.init("--test")

test = code_aster.TestCase()

print("starting...", flush=True)
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    handled = False
    try:
        raise ValueError("first error")
    except ValueError:
        handled = True
        print("don't worry!")
    test.assertTrue(handled)
    raise ValueError("second error")
    MPI.COMM_WORLD.send(None, dest=1, tag=42)
else:
    print("waiting for #0...", flush=True)
    comm.recv(source=0, tag=42)

test.printSummary()

code_aster.close()
