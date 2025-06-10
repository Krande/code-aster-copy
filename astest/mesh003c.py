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

rank = MPI.ASTER_COMM_WORLD.Get_rank()
print("Nb procs", MPI.ASTER_COMM_WORLD.Get_size())
print("Rank", MPI.ASTER_COMM_WORLD.Get_rank())

from code_aster.MedUtils import splitMeshAndFieldsFromMedFile

ret = splitMeshAndFieldsFromMedFile("fort.20", deterministic=True)
pMesh = ret[0]


# Test full mesh
test.assertEqual(pMesh.getDimension(), 3)


# Test ConnectionMesh - The full mesh
print("cMesh1", flush=True)
cMesh1 = CA.ConnectionMesh(pMesh, ["CUBE"], [])
test.assertEqual(cMesh1.getParallelMesh().getName(), pMesh.getName())
test.assertEqual(cMesh1.getDimension(), 3)
test.assertEqual(cMesh1.getNumberOfNodes(), 197)
test.assertEqual(cMesh1.getNumberOfCells(), 973)
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfNodes()), ["CUBE"])
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfCells()), [])
test.assertTrue(cMesh1.hasGroupOfNodes("CUBE"))
test.assertFalse(cMesh1.hasGroupOfNodes("UCBE"))
test.assertFalse(cMesh1.hasGroupOfCells("CUBE"))
test.assertSequenceEqual(sorted(cMesh1.getCells("CUBE")), [])


# Test ConnectionMesh - The full mesh
print("cMesh2", flush=True)
cMesh2 = CA.ConnectionMesh(pMesh, [], ["CUBE"])
test.assertEqual(cMesh2.getDimension(), 3)
test.assertEqual(cMesh2.getNumberOfNodes(), 197)
test.assertEqual(cMesh2.getNumberOfCells(), 973)
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfCells()), ["CUBE"])
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfNodes()), [])
test.assertTrue(cMesh2.hasGroupOfCells("CUBE"))
test.assertFalse(cMesh2.hasGroupOfCells("UCBE"))
test.assertFalse(cMesh2.hasGroupOfNodes("CUBE"))


# Test ConnectionMesh - a part mesh
print("cMesh3", flush=True)
cMesh3 = CA.ConnectionMesh(pMesh, ["N0", "N2"], [])
test.assertEqual(cMesh3.getDimension(), 3)
test.assertEqual(cMesh3.getNumberOfNodes(), 11)
test.assertEqual(cMesh3.getNumberOfCells(), 20)
test.assertSequenceEqual(sorted(cMesh3.getGroupsOfNodes()), ["N0", "N2"])
test.assertSequenceEqual(
    sorted(cMesh3.getCells()),
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
)
test.assertSequenceEqual(sorted(cMesh3.getCells("CUBE")), [])
test.assertSequenceEqual(
    cMesh3.getNodesGlobalNumbering(), [3, 4, 16, 23, 27, 28, 44, 51, 77, 103, 135]
)
if rank == 0:
    test.assertSequenceEqual(
        sorted(cMesh3.getNodesLocalNumbering()), [-1, -1, -1, -1, -1, 4, 10, 11, 25, 33, 55]
    )
elif rank == 1:
    test.assertSequenceEqual(
        sorted(cMesh3.getNodesLocalNumbering()), [-1, -1, -1, -1, -1, -1, 31, 39, 46, 51, 68]
    )
elif rank == 2:
    test.assertSequenceEqual(
        sorted(cMesh3.getNodesLocalNumbering()), [-1, -1, -1, -1, -1, -1, -1, -1, -1, 57, 79]
    )
else:
    assert False

# Test ConnectionMesh - a part mesh
print("cMesh4", flush=True)
cMesh4 = CA.ConnectionMesh(pMesh, ["ALL_NO"], [])
test.assertEqual(cMesh4.getDimension(), 3)
test.assertEqual(cMesh4.getNumberOfNodes(), 41)
test.assertEqual(cMesh4.getNumberOfCells(), 74)
test.assertSequenceEqual(sorted(cMesh4.getGroupsOfNodes()), ["ALL_NO"])


# Test ConnectionMesh - a part mesh
print("cMesh5", flush=True)
cMesh5 = CA.ConnectionMesh(pMesh, [], ["ALL_SEG"])
test.assertEqual(cMesh5.getDimension(), 3)
test.assertEqual(cMesh5.getNumberOfNodes(), 139)
test.assertEqual(cMesh5.getNumberOfCells(), 422)
test.assertSequenceEqual(sorted(cMesh5.getGroupsOfCells()), ["ALL_SEG"])


# Test ConnectionMesh - a part mesh
print("cMesh6", flush=True)
cMesh6 = CA.ConnectionMesh(pMesh, [], ["S1", "S2"])
test.assertEqual(cMesh6.getDimension(), 3)
test.assertEqual(cMesh6.getNumberOfNodes(), 41)
test.assertEqual(cMesh6.getNumberOfCells(), 103)
test.assertSequenceEqual(sorted(cMesh6.getGroupsOfCells()), ["S1", "S2"])
test.assertSequenceEqual(sorted(cMesh6.getGroupsOfNodes()), [])


# Test ConnectionMesh - a part mesh
print("cMesh7", flush=True)
cMesh7 = CA.ConnectionMesh(pMesh, ["N1", "N4"], ["S1", "S2"])
test.assertEqual(cMesh7.getDimension(), 3)
test.assertEqual(cMesh7.getNumberOfNodes(), 52)
test.assertEqual(cMesh7.getNumberOfCells(), 123)
test.assertSequenceEqual(sorted(cMesh6.getGroupsOfCells()), ["S1", "S2"])
test.assertSequenceEqual(sorted(cMesh7.getGroupsOfNodes()), ["N1", "N4"])


# Test ConnectionMesh - a part mesh
print("cMesh8", flush=True)
cMesh8 = CA.ConnectionMesh(pMesh, [], ["OUEST", "NORD"])
test.assertEqual(cMesh8.getDimension(), 3)
test.assertEqual(cMesh8.getNumberOfNodes(), 116)
test.assertEqual(cMesh8.getNumberOfCells(), 460)
test.assertSequenceEqual(sorted(cMesh8.getGroupsOfCells()), ["NORD", "OUEST"])


# Test ConnectionMesh - a part mesh
print("cMesh9", flush=True)
cMesh9 = CA.ConnectionMesh(pMesh, [], ["OUEST", "NORD", "S1"])
test.assertEqual(cMesh9.getDimension(), 3)
test.assertEqual(cMesh9.getNumberOfNodes(), 129)
test.assertEqual(cMesh9.getNumberOfCells(), 497)
test.assertSequenceEqual(sorted(cMesh9.getGroupsOfCells()), ["NORD", "OUEST", "S1"])


# Test ConnectionMesh - a part mesh
print("cMesh10", flush=True)
cMesh10 = CA.ConnectionMesh(pMesh, ["N1"], ["OUEST"])
test.assertEqual(cMesh10.getDimension(), 3)
test.assertEqual(cMesh10.getNumberOfNodes(), 75)
test.assertEqual(cMesh10.getNumberOfCells(), 270)
test.assertSequenceEqual(sorted(cMesh10.getGroupsOfCells()), ["OUEST"])
test.assertSequenceEqual(sorted(cMesh10.getGroupsOfNodes()), ["N1"])


# Test ConnectionMesh - The full mesh
print("cMesh11", flush=True)
cMesh11 = CA.ConnectionMesh(pMesh, ["CUBE"], ["CUBE"])
test.assertEqual(cMesh11.getDimension(), 3)
test.assertEqual(cMesh11.getNumberOfNodes(), 197)
test.assertEqual(cMesh11.getNumberOfCells(), 973)
test.assertSequenceEqual(sorted(cMesh11.getGroupsOfNodes()), ["CUBE"])
test.assertSequenceEqual(sorted(cMesh11.getGroupsOfCells()), ["CUBE"])
test.assertTrue(cMesh11.hasGroupOfCells("CUBE"))
test.assertTrue(cMesh11.hasGroupOfNodes("CUBE"))
test.assertFalse(cMesh11.hasGroupOfCells("UCBE"))
test.assertFalse(cMesh11.hasGroupOfNodes("UCBE"))


# Test ConnectionMesh - The full mesh
print("cMesh12", flush=True)
cMesh12 = CA.ConnectionMesh(pMesh, ["ALL_NO", "N0", "N1", "N2"], ["ALL_SEG", "S1", "S2"])
test.assertEqual(cMesh12.getDimension(), 3)
test.assertEqual(cMesh12.getNumberOfNodes(), 139)
test.assertEqual(cMesh12.getNumberOfCells(), 422)
test.assertSequenceEqual(sorted(cMesh12.getGroupsOfNodes()), ["ALL_NO", "N0", "N1", "N2"])
test.assertSequenceEqual(sorted(cMesh12.getGroupsOfCells()), ["ALL_SEG", "S1", "S2"])


# Test ConnectionMesh
print("cMesh13", flush=True)
cMesh13 = CA.ConnectionMesh(
    pMesh, ["N0", "N1", "N2", "N3"], ["ALL_SEG", "BAS", "HAUT", "SUD", "XXX"]
)
test.assertEqual(cMesh13.getParallelMesh().getName(), pMesh.getName())
test.assertEqual(cMesh13.getDimension(), 3)
test.assertEqual(cMesh13.getNumberOfNodes(), 174)
test.assertEqual(cMesh13.getNumberOfCells(), 710)
test.assertSequenceEqual(sorted(cMesh13.getGroupsOfNodes()), ["N0", "N1", "N2"])
test.assertSequenceEqual(sorted(cMesh13.getGroupsOfCells()), ["ALL_SEG", "BAS", "HAUT", "SUD"])

# Test model
cModel = AFFE_MODELE(
    MAILLAGE=cMesh13,
    AFFE=(
        _F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"),
        _F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="DIS_T"),
    ),
    VERI_JACOBIEN="NON",
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)

test.assertEqual(cMesh13.getName(), cModel.getMesh().getName())

test.printSummary()


FIN()
