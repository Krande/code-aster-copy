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

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

from code_aster.MedUtils import splitMeshAndFieldsFromMedFile

ret = splitMeshAndFieldsFromMedFile("fort.20", deterministic=True)
pMesh = ret[0]

# Test node gathering while partitioning mesh
ret2 = splitMeshAndFieldsFromMedFile("fort.20", deterministic=True, nodeGrpToGather=[["ToGather"]])
pMesh2 = ret2[0]
nodeList = pMesh2.getNodes("ToGather")
nodeOwner = pMesh2.getNodesOwner()
if len(nodeList) != 0:
    masterNodeOwner = nodeOwner[nodeList[0]]
    # Test if there is only one proc number for all node group
    for nodeId in nodeList:
        test.assertEqual(nodeOwner[nodeId], masterNodeOwner)
else:
    test.assertTrue(True)

model = AFFE_MODELE(
    MAILLAGE=pMesh, AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI")
)

rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()

if nbproc == 2:
    nbNodes = [118, 112]
    nbCells = [56, 104]
elif nbproc == 3:
    nbNodes = [72, 103, 80]
    nbCells = [56, 48, 67]
elif nbproc == 4:
    nbNodes = [61, 85, 60, 76]
    nbCells = [31, 31, 37, 78]

test.assertEqual(pMesh.getDimension(), 2)
test.assertEqual(pMesh.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(pMesh.getNumberOfCells(), nbCells[rank])
test.assertTrue(pMesh.isParallel())

test.printSummary()

FIN()
