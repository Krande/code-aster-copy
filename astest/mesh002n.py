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
from code_aster import MPI
from code_aster.Utilities import config

code_aster.init("--test")

test = code_aster.TestCase()

pMesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH")

model = AFFE_MODELE(MAILLAGE=pMesh, AFFE=_F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"))

rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()

if nbproc == 1:
    nbNodes = [778]
    nbCells = [647]
elif nbproc == 2:
    nbNodes = [454, 498]
    nbCells = [350, 388]
elif nbproc == 3:
    nbNodes = [340, 446, 344]
    nbCells = [263, 305, 266]
elif nbproc == 4:
    nbNodes = [307, 314, 328, 364]
    nbCells = [213, 222, 236, 253]

test.assertEqual(pMesh.getDimension(), 3)
test.assertEqual(pMesh.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(pMesh.getNumberOfCells(), nbCells[rank])
test.assertTrue(pMesh.isParallel() or not config.get("parallel"))

test.printSummary()


FIN()
