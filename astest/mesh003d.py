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
cMesh1 = ConnectionMesh(pMesh, [], ["COMPLET"])
test.assertEqual(cMesh1.getParallelMesh().getName(), pMesh.getName())
test.assertEqual(cMesh1.getDimension(), 3)
test.assertEqual(cMesh1.getNumberOfNodes(), 50193)
test.assertEqual(cMesh1.getNumberOfCells(), 25047)
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfNodes()), [])
test.assertSequenceEqual(sorted(cMesh1.getGroupsOfCells()), ["COMPLET"])
test.assertFalse( cMesh1.hasGroupOfNodes("COMPLET") )
test.assertFalse( cMesh1.hasGroupOfCells("AFCE") )
test.assertTrue( cMesh1.hasGroupOfCells("COMPLET") )


# Test ConnectionMesh - The full mesh
print("cMesh2", flush=True)
cMesh2 = ConnectionMesh(pMesh, [], ["PIPEFISS"])
test.assertEqual(cMesh2.getDimension(), 3)
test.assertEqual(cMesh2.getNumberOfNodes(), 23731)
test.assertEqual(cMesh2.getNumberOfCells(), 9938)
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfCells()), ["PIPEFISS"])
test.assertSequenceEqual(sorted(cMesh2.getGroupsOfNodes()), [])
test.assertTrue( cMesh2.hasGroupOfCells("PIPEFISS") )
test.assertFalse( cMesh2.hasGroupOfCells("AFCE") )
test.assertFalse( cMesh2.hasGroupOfNodes("FACE") )

# Test ConnectionMesh - a part mesh
print("cMesh3", flush=True)
cMesh3 = ConnectionMesh(pMesh, [], ["facePeau0", "facePeau1", "haut"])
test.assertEqual(cMesh3.getDimension(), 3)
test.assertEqual(cMesh3.getNumberOfNodes(), 7101)
test.assertEqual(cMesh3.getNumberOfCells(), 4024)
test.assertSequenceEqual(sorted(cMesh3.getGroupsOfCells()), ["facePeau0", "facePeau1", "haut"])


# Test ConnectionMesh - a part mesh
print("cMesh4", flush=True)
cMesh4 = ConnectionMesh(pMesh, ["FONDFISS"], [])
test.assertEqual(cMesh4.getDimension(), 3)
test.assertEqual(cMesh4.getNumberOfNodes(), 3158)
test.assertEqual(cMesh4.getNumberOfCells(), 869)
test.assertSequenceEqual(sorted(cMesh4.getGroupsOfNodes()), ["FONDFISS"])


# Test ConnectionMesh - a part mesh
print("cMesh5", flush=True)
cMesh5 = ConnectionMesh(pMesh, [], ["bords"])
test.assertEqual(cMesh5.getDimension(), 3)
test.assertEqual(cMesh5.getNumberOfNodes(), 426)
test.assertEqual(cMesh5.getNumberOfCells(), 192)
test.assertSequenceEqual(sorted(cMesh5.getGroupsOfCells()), ["bords"])
test.assertSequenceEqual(sorted(cMesh5.getCells("bords")), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])


# Test ConnectionMesh - Isolated node
print("cMesh6", flush=True)
cMesh6 = ConnectionMesh(pMesh, ["nfondfis"], ["bas", "haut", "bords", "affVols"])
test.assertEqual(cMesh6.getDimension(), 3)
test.assertEqual(cMesh6.getNumberOfNodes(), 8928)
test.assertEqual(cMesh6.getNumberOfCells(), 3486)
test.assertSequenceEqual(sorted(cMesh6.getGroupsOfNodes()), ["nfondfis"])
test.assertSequenceEqual(sorted(cMesh6.getGroupsOfCells()), sorted(["bas", "haut", "bords", "affVols"]))
test.assertEqual(sum(list(cMesh6.getGlobalNumbering())), 192089986)
test.assertEqual(sum(list(cMesh6.getLocalNumbering())), 52960885)


# Test model
cModel = AFFE_MODELE(MAILLAGE=cMesh6,
                    AFFE=(_F(TOUT='OUI', PHENOMENE='MECANIQUE',
                                         MODELISATION="D_PLAN",),
                          _F(TOUT='OUI', PHENOMENE='MECANIQUE',
                                         MODELISATION='DIS_T',),
                        ),
                    VERI_JACOBIEN='NON',
                    DISTRIBUTION=_F(METHODE='CENTRALISE',),)

test.assertEqual(cMesh6.getName(), cModel.getMesh().getName())

test.printSummary()


FIN()
