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
import numpy as np

from code_aster.Commands import DEFI_GROUP
from code_aster import MPI


code_aster.init("--test")

# check ParallelMesh object API
test = code_aster.TestCase()

# MPI test
rank = MPI.COMM_WORLD.Get_rank()
nbproc = MPI.COMM_WORLD.Get_size()

test.assertEqual(nbproc, 3)

# from MED format (only this one a ParallelMesh)
mesh = code_aster.ParallelMesh()
mesh.readMedFile("mesh001b/%d.med" % rank, True)


mesh = DEFI_GROUP(reuse=mesh, MAILLAGE=mesh,
                  CREA_GROUP_NO=(_F(NOM='GN'+str(rank),  GROUP_NO='EXT_0'),),
                  CREA_GROUP_MA=(_F(NOM='GC'+str(rank),  GROUP_MA='Cable0'),))

test.assertTrue(mesh.isParallel())
test.assertEqual(mesh.getDimension(), 3)

nbNodes = [89, 90, 109]
nbCells = [59, 54, 75]

test.assertEqual(mesh.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(mesh.getNumberOfCells(), nbCells[rank])

# test groups
globalGroupsOfCells = ["Cable0", "Cable1", "Cable2", "Cable3", "Cable4", "Cable5", "Cable6",
                       "Cable7", "Cable8", "Press", "Encast", "Beton", "GC0", "GC1", "GC2"]
groupsOfCells = [["Cable0", "Cable1", "Cable2", "Cable3", "Cable4", "Cable5", "Cable6",
                  "Cable7", "Cable8", "Press", "Encast", "Beton", "GC0"],
                 ["Cable0", "Cable1", "Cable2", "Cable3", "Cable4", "Cable5", "Cable6",
                  "Cable7", "Cable8", "Press", "Beton", "GC1"],
                 ["Cable0", "Cable1", "Cable2", "Cable3", "Cable4", "Cable5", "Cable6",
                  "Cable7", "Cable8", "Press", "Encast", "Beton", "GC2"]]

test.assertSequenceEqual(sorted(mesh.getGroupsOfCells()),
                         sorted(globalGroupsOfCells))
test.assertSequenceEqual(
    sorted(mesh.getGroupsOfCells(False)), sorted(globalGroupsOfCells))
test.assertSequenceEqual(
    sorted(mesh.getGroupsOfCells(True)), sorted(groupsOfCells[rank]))

test.assertTrue(mesh.hasGroupOfCells("Beton"))
test.assertTrue(mesh.hasGroupOfCells("Beton", False))
test.assertTrue(mesh.hasGroupOfCells("Beton", True))
test.assertTrue(mesh.hasGroupOfCells("GC1", False))
test.assertTrue(mesh.hasGroupOfCells("GC"+str(rank), True))
test.assertFalse(mesh.hasGroupOfCells("GC4", True))
test.assertFalse(mesh.hasGroupOfCells("GC4", False))

globalGroupsOfNodes = ["EXT_0", "EXT_1", "EXT_2", "GN0", "GN1", "GN2"]
groupsOfNodes = [["EXT_0", "EXT_1", "EXT_2", "GN0"], ["EXT_0", "EXT_1", "EXT_2", "GN1"],
                 ["EXT_0", "EXT_1", "EXT_2", "GN2"]]

test.assertSequenceEqual(sorted(mesh.getGroupsOfNodes()),
                         sorted(globalGroupsOfNodes))
test.assertSequenceEqual(
    sorted(mesh.getGroupsOfNodes(False)), sorted(globalGroupsOfNodes))
test.assertSequenceEqual(
    sorted(mesh.getGroupsOfNodes(True)), sorted(groupsOfNodes[rank]))

test.assertTrue(mesh.hasGroupOfNodes("EXT_0"))
test.assertTrue(mesh.hasGroupOfNodes("EXT_0", False))
test.assertTrue(mesh.hasGroupOfNodes("EXT_0", True))
test.assertTrue(mesh.hasGroupOfNodes("GN1", False))
test.assertTrue(mesh.hasGroupOfNodes("GN"+str(rank), True))
test.assertFalse(mesh.hasGroupOfNodes("GN4", True))
test.assertFalse(mesh.hasGroupOfNodes("GN4", False))

# Link between local and global numbering
globalNodesNum = mesh.getNodes(localNumbering=False)
nodesGlobFirst = [4, 44, 0]
test.assertEqual(globalNodesNum[0], nodesGlobFirst[rank])
nodesGlobLast = [97, 97, 95]
test.assertEqual(globalNodesNum[-1], nodesGlobLast[rank])

# Owner of Nodes
NodesRank = mesh.getNodesRank()
test.assertEqual(NodesRank[1], rank)

# Node 92 (index is 91) is shared by all meshes (owner is 1)
node92 = [85, 84, 106]
test.assertEqual(globalNodesNum[node92[rank]-1], 92-1)
test.assertEqual(NodesRank[node92[rank]-1], 1)

# Cells rank
cellsRankRef = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0, 1, 1, 0, 0, 1,
                                                             1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0,
                                                             1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0],
                [2, 2, 0, 2, 1, 1, 2, 1, 2, 2, 0, 2, 1, 1, 2, 1, 2, 2, 0, 2, 1,
                 0, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2,
                 1, 2, 2, 2, 0, 2, 2, 1, 1, 2, 2, 1, 2, 2, 2, 0, 2, 2, 1,
                 1, 2, 2, 1, 2, 2, 2, 0, 2, 2, 1, 0, 2, 2, 1]]
test.assertEqual(len(mesh.getCellsRank()), mesh.getNumberOfCells())
test.assertSequenceEqual(mesh.getCellsRank(), cellsRankRef[rank])


def inter(list1, list2):
    return list(set(list1).intersection(list2))


test.assertEqual(mesh.getNumberOfNodes(), len(mesh.getNodes()))
test.assertEqual(mesh.getNumberOfNodes(), len(
    mesh.getNodes(localNumbering=True)))
test.assertEqual(mesh.getNumberOfNodes(), len(
    mesh.getNodes(localNumbering=False)))

test.assertEqual(mesh.getNumberOfCells(), len(mesh.getCells()))

test.assertSequenceEqual(mesh.getNodes(), range(1, mesh.getNumberOfNodes()+1))
test.assertSequenceEqual(mesh.getNodes(
    localNumbering=True), range(1, mesh.getNumberOfNodes()+1))
test.assertSequenceEqual(mesh.getCells(), range(1, mesh.getNumberOfCells()+1))

allNodes = []
innerNodes = []
outerNodes = []

for i in range(mesh.getNumberOfNodes()):
    allNodes.append(i+1)
    if NodesRank[i] == rank:
        innerNodes.append(i+1)
    else:
        outerNodes.append(i+1)

test.assertTrue(len(inter(innerNodes, outerNodes)) == 0)


test.assertSequenceEqual(sorted(mesh.getNodes()), sorted(allNodes))
test.assertSequenceEqual(
    sorted(mesh.getNodes(localNumbering=True)), sorted(allNodes))
test.assertSequenceEqual(sorted(mesh.getNodes(localNumbering=False)),
                         sorted([globalNodesNum[i-1] for i in allNodes]))
test.assertSequenceEqual(sorted(mesh.getNodes(
    localNumbering=True, same_rank=False)), sorted(outerNodes))
test.assertSequenceEqual(sorted(mesh.getNodes(localNumbering=False, same_rank=False)),
                         sorted([globalNodesNum[i-1] for i in outerNodes]))
test.assertSequenceEqual(sorted(mesh.getNodes(
    localNumbering=True, same_rank=True)), sorted(innerNodes))
test.assertSequenceEqual(sorted(mesh.getNodes(localNumbering=False, same_rank=True)),
                         sorted([globalNodesNum[i-1] for i in innerNodes]))

test.assertSequenceEqual(sorted(mesh.getNodes("Beton")),
                         sorted(mesh.getNodes("Beton", True, None)))
test.assertSequenceEqual(sorted(mesh.getNodes("Beton")),
                         sorted(inter(mesh.getNodes("Beton"), mesh.getNodes())))
test.assertSequenceEqual(sorted(mesh.getNodes("Beton")),
                         sorted(inter(mesh.getNodes("Beton", True), mesh.getNodes(localNumbering=True))))
test.assertSequenceEqual(sorted(mesh.getNodes("Beton")),
                         sorted(inter(mesh.getNodes("Beton"), mesh.getNodes(localNumbering=True))))
test.assertSequenceEqual(sorted(mesh.getNodes("Beton")),
                         sorted(inter(mesh.getNodes("Beton", True), mesh.getNodes())))
test.assertSequenceEqual(sorted(mesh.getNodes("Beton", False)),
                         sorted(inter(mesh.getNodes("Beton", False), mesh.getNodes(localNumbering=False))))
test.assertSequenceEqual(sorted(mesh.getNodes("Beton", False)),
                         sorted(inter(mesh.getNodes("Beton", False), mesh.getNodes(localNumbering=False))))
test.assertSequenceEqual(sorted(mesh.getNodes("Beton", False)),
                         sorted(inter(mesh.getNodes("Beton", False), allNodes)))

# rank-owned nodes in local numbering
nodeLoc = [[59, 69, 70], [], [77, 86, 87]]
test.assertSequenceEqual(mesh.getNodesFromCells("Cable0", True, True),
                         nodeLoc[rank])
test.assertSequenceEqual(mesh.getNodesFromCells("Cable0", False, True),
                         sorted([globalNodesNum[i-1] for i in nodeLoc[rank]]))

# nodes in local numbering
nodeLocWithGhosts = [[59, 68, 69, 70], [61, 62], [77, 86, 87, 88]]
test.assertSequenceEqual(mesh.getNodesFromCells("Cable0", True, None),
                         nodeLocWithGhosts[rank])
test.assertSequenceEqual(mesh.getNodesFromCells("Cable0", False, None),
                         sorted([globalNodesNum[i-1] for i in nodeLocWithGhosts[rank]]))


new_mesh = mesh.refine(2)

# test mesh builder
test.assertEqual(code_aster.ParallelMesh.buildSquare(refine=4).getNumberOfNodes(),
                 [116, 123, 126][rank])
test.assertEqual(code_aster.ParallelMesh.buildCube(refine=4).getNumberOfNodes(),
                 [1929, 2079, 2189][rank])
test.assertEqual(code_aster.ParallelMesh.buildDisk(refine=5).getNumberOfNodes(),
                 [4114, 4421, 4404][rank])
test.assertEqual(code_aster.ParallelMesh.buildCylinder(refine=3).getNumberOfNodes(),
                 [5246, 5242, 5625][rank])

test.printSummary()

code_aster.close()
