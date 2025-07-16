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


CA.init("--test")

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
cMesh1 = CA.ConnectionMesh(pMesh, [], ["VTOT"])
test.assertEqual(cMesh1.getParallelMesh().getName(), pMesh.getName())
test.assertEqual(cMesh1.getDimension(), 3)
test.assertEqual(cMesh1.getNumberOfNodes(), 17331)
test.assertEqual(cMesh1.getNumberOfCells(), 8798)
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfNodes()), [])
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfCells()), ["VTOT"])
test.assertFalse(cMesh1.hasGroupOfNodes("AFCE"))
test.assertTrue(cMesh1.hasGroupOfCells("VTOT"))
test.assertEqual(sum(list(cMesh1.getCells())), 38698003)
test.assertEqual(sum(list(cMesh1.getNodesGlobalNumbering())), 150173315)
if rank == 0:
    test.assertEqual(sum(cMesh1.getNodesLocalNumbering()), 16846196)
elif rank == 1:
    test.assertEqual(sum(cMesh1.getNodesLocalNumbering()), 15671069)
elif rank == 2:
    test.assertEqual(sum(cMesh1.getNodesLocalNumbering()), 17550294)
else:
    assert False


# Test ConnectionMesh - The full mesh
print("cMesh2", flush=True)
cMesh2 = CA.ConnectionMesh(pMesh, ["A9", "D7", "B3"], ["CD9", "AB5", "AB1", "ETE4"])
test.assertEqual(cMesh2.getDimension(), 3)
test.assertEqual(cMesh2.getNumberOfNodes(), 826)
test.assertEqual(cMesh2.getNumberOfCells(), 1295)
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfCells()), sorted(["CD9", "AB5", "AB1", "ETE4"]))
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfNodes()), sorted(["A9", "D7", "B3"]))
test.assertTrue(cMesh2.hasGroupOfCells("AB1"))
test.assertFalse(cMesh2.hasGroupOfCells("AFCE"))
test.assertFalse(cMesh2.hasGroupOfNodes("FACE"))
test.assertEqual(sum(list(cMesh2.getCells())), 837865)
print(len(cMesh2.getCells("AB1")))
test.assertEqual(sum(list(cMesh2.getCells("AB1"))), 23630)
test.assertEqual(sum(list(cMesh2.getNodesGlobalNumbering())), 5193300)
if rank == 0:
    test.assertEqual(sum(cMesh2.getNodesLocalNumbering()), 477000)
elif rank == 1:
    test.assertEqual(sum(cMesh2.getNodesLocalNumbering()), 21026)
elif rank == 2:
    test.assertEqual(sum(cMesh2.getNodesLocalNumbering()), 1352010)
else:
    assert False

test.printSummary()


FIN()
