# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
from code_aster.Objects import ConnectionMesh
from code_aster.Commands import *
from code_aster import MPI


code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.COMM_WORLD.Get_rank()
print("Nb procs", MPI.COMM_WORLD.Get_size())
print("Rank", MPI.COMM_WORLD.Get_rank())

pMesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH")


# Test full mesh
test.assertEqual(pMesh.getDimension(), 3)


# Test ConnectionMesh - The full mesh
print("cMesh1", flush=True)
cMesh1 = ConnectionMesh(pMesh, [], ["VTOT"])
test.assertEqual(cMesh1.getParallelMesh().getName(), pMesh.getName())
test.assertEqual(cMesh1.getDimension(), 3)
test.assertEqual(cMesh1.getNumberOfNodes(), 17331)
test.assertEqual(cMesh1.getNumberOfCells(), 8798)
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfNodes()), [])
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfCells()), ["VTOT"])
test.assertFalse( cMesh1.hasGroupOfNodes("AFCE") )
test.assertTrue( cMesh1.hasGroupOfCells("VTOT") )
test.assertEqual(sum(list(cMesh1.getCells())), 38706801)
test.assertEqual(sum(list(cMesh1.getGlobalNumbering())), 150173315)
test.assertEqual(sum(list(cMesh1.getLocalNumbering())), 50096296)


# Test ConnectionMesh - The full mesh
print("cMesh2", flush=True)
cMesh2 = ConnectionMesh(pMesh, ["A9", "D7", "B3"], ["CD9", "AB5", "AB1", "ETE4"])
test.assertEqual(cMesh2.getDimension(), 3)
test.assertEqual(cMesh2.getNumberOfNodes(), 826)
test.assertEqual(cMesh2.getNumberOfCells(), 1295)
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfCells()), sorted(["CD9", "AB5", "AB1", "ETE4"]))
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfNodes()), sorted(["A9", "D7", "B3"]))
test.assertTrue( cMesh2.hasGroupOfCells("AB1") )
test.assertFalse( cMesh2.hasGroupOfCells("AFCE") )
test.assertFalse( cMesh2.hasGroupOfNodes("FACE") )
test.assertEqual(sum(list(cMesh2.getCells())), 839160)
test.assertEqual(sum(list(cMesh2.getCells("AB1"))), 23650)
test.assertEqual(sum(list(cMesh2.getGlobalNumbering())), 5193300)
test.assertEqual(sum(list(cMesh2.getLocalNumbering())), 1851292)

test.printSummary()


FIN()
