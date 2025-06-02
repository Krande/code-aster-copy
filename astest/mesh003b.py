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
test.assertEqual(pMesh.getDimension(), 2)


# Test ConnectionMesh - The full mesh
print("cMesh1", flush=True)
cMesh1 = CA.ConnectionMesh(pMesh, ["FACE"], [])
test.assertEqual(cMesh1.getParallelMesh().getName(), pMesh.getName())
test.assertEqual(cMesh1.getDimension(), 2)
test.assertEqual(cMesh1.getNumberOfNodes(), 341)
test.assertEqual(cMesh1.getNumberOfCells(), 115)
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfNodes()), ["FACE"])
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfCells()), [])
test.assertTrue(cMesh1.hasGroupOfNodes("FACE"))
test.assertFalse(cMesh1.hasGroupOfNodes("AFCE"))
test.assertFalse(cMesh1.hasGroupOfCells("FACE"))
test.assertSequenceEqual(sorted(cMesh1.getCells()), [i for i in range(115)])
test.assertSequenceEqual(sorted(cMesh1.getCells("FACE")), [])
test.assertSequenceEqual(sorted(cMesh1.getNodesGlobalNumbering()), [i for i in range(0, 341)])

# Test ConnectionMesh - The full mesh
print("cMesh2", flush=True)
cMesh2 = CA.ConnectionMesh(pMesh, [], ["FACE"])
test.assertEqual(cMesh2.getDimension(), 2)
test.assertEqual(cMesh2.getNumberOfNodes(), 341)
test.assertEqual(cMesh2.getNumberOfCells(), 115)
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfCells()), ["FACE"])
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfNodes()), [])
test.assertTrue(cMesh2.hasGroupOfCells("FACE"))
test.assertFalse(cMesh2.hasGroupOfCells("AFCE"))
test.assertFalse(cMesh2.hasGroupOfNodes("FACE"))
test.assertSequenceEqual(sorted(cMesh2.getCells()), [i for i in range(115)])
test.assertEqual(len(cMesh2.getCells("FACE")), 75)
test.assertSequenceEqual(sorted(cMesh2.getNodesGlobalNumbering()), [i for i in range(0, 341)])
test.assertEqual(len(cMesh2.getNodesLocalNumbering()), 341)


# Test ConnectionMesh - a part mesh
print("cMesh3", flush=True)
cMesh3 = CA.ConnectionMesh(pMesh, ["N0", "N2"], [])
test.assertEqual(cMesh3.getDimension(), 2)
test.assertEqual(cMesh3.getNumberOfNodes(), 18)
test.assertEqual(cMesh3.getNumberOfCells(), 6)
test.assertSequenceEqual(sorted(cMesh3.getGroupsOfNodes()), ["N0", "N2"])
test.assertSequenceEqual(sorted(cMesh3.getCells()), [0, 1, 2, 3, 4, 5])
test.assertSequenceEqual(sorted(cMesh3.getCells("FACE")), [])
test.assertSequenceEqual(
    sorted(cMesh3.getNodesGlobalNumbering()),
    [1, 2, 7, 21, 22, 26, 40, 95, 100, 115, 116, 121, 136, 137, 138, 310, 338, 340],
)
if rank == 0:
    test.assertSequenceEqual(
        sorted(cMesh3.getNodesLocalNumbering()),
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, 2, 6, 12, 17, 41, 47, 52, 53, 54],
    )
elif rank == 1:
    test.assertSequenceEqual(
        sorted(cMesh3.getNodesLocalNumbering()),
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    )
elif rank == 2:
    test.assertSequenceEqual(
        sorted(cMesh3.getNodesLocalNumbering()),
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, 35, 46, 48, 49, 56, 57, 89, 95, 96],
    )
else:
    assert False


# Test ConnectionMesh - a part mesh
print("cMesh4", flush=True)
cMesh4 = CA.ConnectionMesh(pMesh, ["N1", "N3"], [])
test.assertEqual(cMesh4.getDimension(), 2)
test.assertEqual(cMesh4.getNumberOfNodes(), 18)
test.assertEqual(cMesh4.getNumberOfCells(), 6)
test.assertSequenceEqual(sorted(cMesh4.getGroupsOfNodes()), ["N1", "N3"])
test.assertSequenceEqual(
    sorted(cMesh4.getNodesGlobalNumbering()),
    [0, 3, 4, 8, 25, 39, 53, 82, 96, 101, 120, 135, 176, 178, 179, 268, 312, 313],
)


# Test ConnectionMesh - a part mesh
print("cMesh5", flush=True)
cMesh5 = CA.ConnectionMesh(pMesh, ["N1", "FACE"], [])
test.assertEqual(cMesh5.getDimension(), 2)
test.assertEqual(cMesh5.getNumberOfNodes(), 341)
test.assertEqual(cMesh5.getNumberOfCells(), 115)
test.assertSequenceEqual(sorted(cMesh5.getGroupsOfNodes()), ["FACE", "N1"])


# Test ConnectionMesh - a part mesh
print("cMesh6", flush=True)
cMesh6 = CA.ConnectionMesh(pMesh, ["N0"], [])
test.assertEqual(cMesh6.getDimension(), 2)
test.assertEqual(cMesh6.getNumberOfNodes(), 9)
test.assertEqual(cMesh6.getNumberOfCells(), 3)
test.assertSequenceEqual(sorted(cMesh6.getGroupsOfNodes()), ["N0"])
test.assertSequenceEqual(
    sorted(cMesh6.getNodesGlobalNumbering()), [1, 7, 26, 40, 100, 121, 136, 137, 138]
)


# Test ConnectionMesh - a part mesh
print("cMesh7", flush=True)
cMesh7 = CA.ConnectionMesh(pMesh, ["N0", "N1", "N2", "N3"], [])
test.assertEqual(cMesh7.getDimension(), 2)
test.assertEqual(cMesh7.getNumberOfNodes(), 36)
test.assertEqual(cMesh7.getNumberOfCells(), 12)
test.assertSequenceEqual(sorted(cMesh7.getGroupsOfNodes()), ["N0", "N1", "N2", "N3"])

# Test ConnectionMesh - a part mesh
print("cMesh8", flush=True)
cMesh8 = CA.ConnectionMesh(pMesh, [], ["GAUCHE"])
test.assertEqual(cMesh8.getDimension(), 2)
test.assertEqual(cMesh8.getNumberOfNodes(), 33)
test.assertEqual(cMesh8.getNumberOfCells(), 12)
test.assertSequenceEqual(sorted(cMesh8.getGroupsOfCells()), ["GAUCHE"])


# Test ConnectionMesh - a part mesh
print("cMesh9", flush=True)
cMesh9 = CA.ConnectionMesh(pMesh, [], ["DROITE"])
test.assertEqual(cMesh9.getDimension(), 2)
test.assertEqual(cMesh9.getNumberOfNodes(), 33)
test.assertEqual(cMesh9.getNumberOfCells(), 12)
test.assertSequenceEqual(sorted(cMesh9.getGroupsOfCells()), ["DROITE"])

# Test ConnectionMesh - a part mesh
print("cMesh10", flush=True)
cMesh10 = CA.ConnectionMesh(pMesh, [], ["GAUCHE", "DROITE"])
test.assertEqual(cMesh10.getDimension(), 2)
test.assertEqual(cMesh10.getNumberOfNodes(), 66)
test.assertEqual(cMesh10.getNumberOfCells(), 24)
test.assertSequenceEqual(sorted(cMesh10.getGroupsOfCells()), ["DROITE", "GAUCHE"])
test.assertEqual(len(cMesh10.getCells()), 24)
test.assertSequenceEqual(sorted(cMesh10.getCells("GAUCHE")), [0, 1, 2, 3, 4])
test.assertSequenceEqual(sorted(cMesh10.getCells("DROITE")), [12, 13, 14, 15, 16])


# Test ConnectionMesh - a part mesh
print("cMesh11", flush=True)
cMesh11 = CA.ConnectionMesh(pMesh, [], ["BAS"])
test.assertEqual(cMesh11.getDimension(), 2)
test.assertEqual(cMesh11.getNumberOfNodes(), 93)
test.assertEqual(cMesh11.getNumberOfCells(), 32)
test.assertSequenceEqual(sorted(cMesh11.getGroupsOfCells()), ["BAS"])


# Test ConnectionMesh - a part mesh
print("cMesh12", flush=True)
cMesh12 = CA.ConnectionMesh(pMesh, [], ["HAUT", "BAS"])
test.assertEqual(cMesh12.getDimension(), 2)
test.assertEqual(cMesh12.getNumberOfNodes(), 186)
test.assertEqual(cMesh12.getNumberOfCells(), 64)
test.assertSequenceEqual(sorted(cMesh12.getGroupsOfCells()), ["BAS", "HAUT"])

# Test ConnectionMesh - The full mesh
print("cMesh13", flush=True)
cMesh13 = CA.ConnectionMesh(pMesh, ["FACE"], ["FACE"])
test.assertEqual(cMesh13.getDimension(), 2)
test.assertEqual(cMesh13.getNumberOfNodes(), 341)
test.assertEqual(cMesh13.getNumberOfCells(), 115)
test.assertSequenceEqual(sorted(cMesh13.getGroupsOfNodes()), ["FACE"])
test.assertSequenceEqual(sorted(cMesh13.getGroupsOfCells()), ["FACE"])
test.assertTrue(cMesh13.hasGroupOfCells("FACE"))
test.assertTrue(cMesh13.hasGroupOfNodes("FACE"))
test.assertFalse(cMesh13.hasGroupOfCells("AFCE"))
test.assertFalse(cMesh13.hasGroupOfNodes("AFCE"))


# Test ConnectionMesh - The full mesh
print("cMesh14", flush=True)
cMesh14 = CA.ConnectionMesh(pMesh, ["N0"], ["DROITE"])
test.assertEqual(cMesh14.getDimension(), 2)
test.assertEqual(cMesh14.getNumberOfNodes(), 42)
test.assertEqual(cMesh14.getNumberOfCells(), 15)
test.assertSequenceEqual(sorted(cMesh14.getGroupsOfNodes()), ["N0"])
test.assertSequenceEqual(sorted(cMesh14.getGroupsOfCells()), ["DROITE"])


# Test ConnectionMesh - The full mesh
print("cMesh15", flush=True)
cMesh15 = CA.ConnectionMesh(pMesh, ["N0", "N1", "N2", "N3"], ["DROITE", "GAUCHE", "BAS", "HAUT"])
test.assertEqual(cMesh15.getParallelMesh().getName(), pMesh.getName())
test.assertEqual(cMesh15.getDimension(), 2)
test.assertEqual(cMesh15.getNumberOfNodes(), 216)
test.assertEqual(cMesh15.getNumberOfCells(), 76)
test.assertSequenceEqual(sorted(cMesh15.getGroupsOfNodes()), ["N0", "N1", "N2", "N3"])
test.assertSequenceEqual(sorted(cMesh15.getGroupsOfCells()), ["BAS", "DROITE", "GAUCHE", "HAUT"])

# Test model
cModel = AFFE_MODELE(
    MAILLAGE=cMesh15,
    AFFE=(
        _F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"),
        _F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="DIS_T"),
    ),
    VERI_JACOBIEN="NON",
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)

test.assertEqual(cMesh15.getName(), cModel.getMesh().getName())

test.printSummary()


FIN()
