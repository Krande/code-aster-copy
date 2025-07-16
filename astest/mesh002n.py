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


from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI
from code_aster.Utilities import config

CA.init("--test")

test = CA.TestCase()

from code_aster.MedUtils import splitMeshAndFieldsFromMedFile

ret = splitMeshAndFieldsFromMedFile("fort.20", deterministic=True)
pMesh = ret[0]

model = AFFE_MODELE(MAILLAGE=pMesh, AFFE=_F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"))

rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()

if nbproc == 1:
    nbNodes = [778]
    nbCells = [647]
elif nbproc == 2:
    nbNodes = [498, 454]
    nbCells = [388, 350]
elif nbproc == 3:
    nbNodes = [355, 418, 344]
    nbCells = [274, 288, 266]
elif nbproc == 4:
    nbNodes = [360, 328, 308, 312]
    nbCells = [249, 236, 216, 217]

test.assertEqual(pMesh.getDimension(), 3)
test.assertEqual(pMesh.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(pMesh.getNumberOfCells(), nbCells[rank])
test.assertTrue(pMesh.isParallel() or not config.get("parallel"))

test.printSummary()


FIN()
