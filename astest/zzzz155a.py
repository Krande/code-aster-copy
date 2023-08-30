# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
from code_aster import MPI
from code_aster.Commands import *

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()


def checkJoints(mesh):
    comm = MPI.COMM_WORLD
    l2G = mesh.getLocalToGlobalMapping()

    j = 0
    for proc in mesh.getOppositeDomains():
        fJ = mesh.getSendJoint(j + 1)
        gFJ = []
        for i in range(int(len(fJ) / 2)):
            gFJ.append(l2G[fJ[2 * i] - 1])

        sJ = mesh.getReceiveJoint(j + 1)
        gSJ = []
        for i in range(int(len(sJ) / 2)):
            gSJ.append(l2G[sJ[2 * i] - 1])

        if proc < rank:
            comm.send(gFJ, dest=proc, tag=j)
            data1 = comm.recv(source=proc, tag=j)
            test.assertEqual(data1 == gFJ, True)
        else:
            data1 = comm.recv(source=proc, tag=j)
            comm.send(gSJ, dest=proc, tag=j)
            test.assertEqual(data1 == gSJ, True)
        j += 1


graph = code_aster.CommGraph()
balancer = code_aster.ObjectBalancer()
a = [i + rank * 10 for i in range(10)]

if rank == 0:
    graph.addCommunication(1)
    graph.addCommunication(3)
    balancer.addElementarySend(1, [0, 2])
    balancer.addElementarySend(3, [1, 3])
    balancer.setElementsToKeep([3])
elif rank == 1:
    graph.addCommunication(0)
    graph.addCommunication(2)
    balancer.addElementarySend(0, [5, 6])
    balancer.addElementarySend(2, [2, 4, 7])
elif rank == 2:
    graph.addCommunication(3)
    balancer.addElementarySend(3, [1, 7])
elif rank == 3:
    graph.addCommunication(1)
    balancer.addElementarySend(1, [8])

balancer.endElementarySendDefinition()
balancer.prepareCommunications()
result = balancer.balanceVectorOverProcesses(a)
print("Result ", result)

if rank == 0:
    test.assertEqual(result[0], 3.0)
    test.assertEqual(result[1], 4.0)
    test.assertEqual(result[2], 5.0)
    test.assertEqual(result[3], 6.0)
    test.assertEqual(result[4], 7.0)
    test.assertEqual(result[5], 8.0)
    test.assertEqual(result[6], 9.0)
    test.assertEqual(result[7], 15.0)
    test.assertEqual(result[8], 16.0)
elif rank == 1:
    test.assertEqual(result[0], 10.0)
    test.assertEqual(result[1], 11.0)
    test.assertEqual(result[2], 13.0)
    test.assertEqual(result[3], 18.0)
    test.assertEqual(result[4], 19.0)
    test.assertEqual(result[5], 0.0)
    test.assertEqual(result[6], 2.0)
    test.assertEqual(result[7], 38.0)
elif rank == 2:
    test.assertEqual(result[0], 20.0)
    test.assertEqual(result[1], 22.0)
    test.assertEqual(result[2], 23.0)
    test.assertEqual(result[3], 24.0)
    test.assertEqual(result[4], 25.0)
    test.assertEqual(result[5], 26.0)
    test.assertEqual(result[6], 28.0)
    test.assertEqual(result[7], 29.0)
    test.assertEqual(result[8], 12.0)
    test.assertEqual(result[9], 14.0)
    test.assertEqual(result[10], 17.0)
elif rank == 3:
    test.assertEqual(result[0], 30.0)
    test.assertEqual(result[1], 31.0)
    test.assertEqual(result[2], 32.0)
    test.assertEqual(result[3], 33.0)
    test.assertEqual(result[4], 34.0)
    test.assertEqual(result[5], 35.0)
    test.assertEqual(result[6], 36.0)
    test.assertEqual(result[7], 37.0)
    test.assertEqual(result[8], 39.0)
    test.assertEqual(result[9], 21.0)
    test.assertEqual(result[10], 27.0)
    test.assertEqual(result[11], 1.0)
    test.assertEqual(result[12], 3.0)

bMesh = code_aster.MeshBalancer()
if rank == 0:
    myMesh = code_aster.Mesh()
    myMesh.readMedFile("fort.20")
    bMesh.buildFromBaseMesh(myMesh)
    outMesh = bMesh.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
elif rank == 1:
    outMesh = bMesh.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
elif rank == 2:
    outMesh = bMesh.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
elif rank == 3:
    outMesh = bMesh.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])

checkMesh = code_aster.Mesh()
checkMesh.readMedFile("fort.20")

mesh2 = code_aster.IncompleteMesh()
mesh2.readMedFile("fort.20")
test.assertEqual(mesh2.debugCheckFromBaseMesh(checkMesh), True)

bMesh = code_aster.MeshBalancer()
bMesh.buildFromBaseMesh(mesh2)
if rank == 0:
    outMesh = bMesh.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
    coords = outMesh.getCoordinates().getValues()
    test.assertEqual(
        coords
        == [
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.0,
            1.0,
            0.0,
            0.0,
            2.0,
            0.0,
            1.0,
            0.0,
            0.0,
            2.0,
            0.0,
            0.0,
            1.0,
            0.0,
            1.0,
            2.0,
            0.0,
            1.0,
            1.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.0,
            1.0,
            2.0,
            0.0,
            2.0,
            2.0,
            0.0,
            1.0,
            2.0,
            1.0,
            2.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            2.0,
            1.0,
            1.0,
        ],
        True,
    )
    connect = outMesh.getConnectivity()
    test.assertEqual(
        connect
        == [
            [2, 1],
            [1, 3],
            [3, 4],
            [2, 5],
            [5, 6],
            [2, 7],
            [7, 8],
            [1, 9],
            [9, 10],
            [2, 1, 3, 5],
            [5, 3, 4, 6],
            [2, 7, 9, 1],
            [7, 8, 10, 9],
            [2, 5, 11, 7],
            [7, 11, 12, 8],
            [5, 6, 13, 11],
            [11, 13, 14, 12],
            [4, 3, 17, 15],
            [15, 17, 18, 16],
            [3, 1, 9, 17],
            [17, 9, 10, 18],
            [1, 2, 5, 3, 9, 7, 11, 17],
            [9, 7, 11, 17, 10, 8, 12, 18],
            [3, 5, 6, 4, 17, 11, 13, 15],
            [17, 11, 13, 15, 18, 12, 14, 16],
        ],
        True,
    )
elif rank == 1:
    outMesh = bMesh.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
    coords = outMesh.getCoordinates().getValues()
    test.assertEqual(
        coords
        == [
            3.0,
            1.0,
            1.0,
            3.0,
            2.0,
            1.0,
            3.0,
            1.0,
            0.0,
            3.0,
            2.0,
            0.0,
            3.0,
            0.0,
            1.0,
            3.0,
            0.0,
            0.0,
            1.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.0,
            1.0,
            2.0,
            0.0,
            2.0,
            2.0,
            0.0,
            1.0,
            2.0,
            1.0,
            2.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            2.0,
            1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            2.0,
            0.0,
            0.0,
            1.0,
            0.0,
            1.0,
            2.0,
            0.0,
            1.0,
        ],
        True,
    )
    connect = outMesh.getConnectivity()
    test.assertEqual(
        connect
        == [
            [6, 5],
            [5, 1],
            [1, 2],
            [6, 3],
            [18, 5],
            [3, 4],
            [15, 16],
            [16, 6],
            [17, 18],
            [15, 16, 18, 17],
            [16, 6, 5, 18],
            [15, 7, 8, 16],
            [5, 6, 3, 1],
            [1, 3, 4, 2],
            [11, 13, 14, 12],
            [12, 14, 1, 2],
            [13, 17, 18, 14],
            [14, 18, 5, 1],
            [16, 8, 3, 6],
            [7, 9, 10, 8],
            [8, 10, 4, 3],
            [18, 16, 8, 14, 5, 6, 3, 1],
            [17, 15, 7, 13, 18, 16, 8, 14],
            [13, 7, 9, 11, 14, 8, 10, 12],
            [14, 8, 10, 12, 1, 3, 4, 2],
        ],
        True,
    )
elif rank == 2:
    outMesh = bMesh.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
    coords = outMesh.getCoordinates().getValues()
    test.assertEqual(
        coords
        == [
            1.0,
            3.0,
            0.0,
            2.0,
            3.0,
            0.0,
            1.0,
            3.0,
            1.0,
            2.0,
            3.0,
            1.0,
            1.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.0,
            1.0,
            2.0,
            0.0,
            2.0,
            2.0,
            0.0,
            1.0,
            2.0,
            1.0,
            2.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            2.0,
            1.0,
            1.0,
            3.0,
            3.0,
            1.0,
            3.0,
            3.0,
            0.0,
            3.0,
            1.0,
            1.0,
            3.0,
            2.0,
            1.0,
            3.0,
            1.0,
            0.0,
            3.0,
            2.0,
            0.0,
        ],
        True,
    )
    connect = outMesh.getConnectivity()
    test.assertEqual(
        connect
        == [
            [17, 18],
            [18, 14],
            [1, 2],
            [2, 14],
            [3, 4],
            [4, 13],
            [15, 16],
            [16, 13],
            [14, 13],
            [5, 7, 8, 6],
            [6, 8, 18, 17],
            [7, 1, 2, 8],
            [8, 2, 14, 18],
            [3, 9, 10, 4],
            [4, 10, 16, 13],
            [9, 11, 12, 10],
            [10, 12, 15, 16],
            [15, 17, 18, 16],
            [16, 18, 14, 13],
            [14, 2, 4, 13],
            [2, 1, 3, 4],
            [11, 5, 7, 9, 12, 6, 8, 10],
            [12, 6, 8, 10, 15, 17, 18, 16],
            [9, 7, 1, 3, 10, 8, 2, 4],
            [10, 8, 2, 4, 16, 18, 14, 13],
        ],
        True,
    )
elif rank == 3:
    outMesh = bMesh.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])
    coords = outMesh.getCoordinates().getValues()
    test.assertEqual(
        coords
        == [
            1.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.0,
            1.0,
            2.0,
            0.0,
            2.0,
            2.0,
            0.0,
            1.0,
            2.0,
            1.0,
            2.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            2.0,
            1.0,
            1.0,
            1.0,
            3.0,
            0.0,
            2.0,
            3.0,
            0.0,
            1.0,
            3.0,
            1.0,
            2.0,
            3.0,
            1.0,
            0.0,
            1.0,
            1.0,
            0.0,
            2.0,
            1.0,
            0.0,
            1.0,
            0.0,
            0.0,
            2.0,
            0.0,
            0.0,
            3.0,
            1.0,
            0.0,
            3.0,
            0.0,
        ],
        True,
    )
    connect = outMesh.getConnectivity()
    test.assertEqual(
        connect
        == [
            [18, 9],
            [9, 10],
            [17, 11],
            [11, 12],
            [16, 18],
            [13, 14],
            [14, 17],
            [18, 17],
            [15, 16],
            [17, 14, 5, 11],
            [11, 5, 6, 12],
            [14, 13, 7, 5],
            [5, 7, 8, 6],
            [15, 16, 3, 1],
            [1, 3, 4, 2],
            [16, 18, 9, 3],
            [3, 9, 10, 4],
            [10, 9, 11, 12],
            [9, 18, 17, 11],
            [15, 13, 14, 16],
            [16, 14, 17, 18],
            [14, 16, 18, 17, 5, 3, 9, 11],
            [5, 3, 9, 11, 6, 4, 10, 12],
            [7, 1, 3, 5, 8, 2, 4, 6],
            [13, 15, 16, 14, 7, 1, 3, 5],
        ],
        True,
    )
checkJoints(outMesh)

part = code_aster.PtScotchPartitioner()
if rank == 0:
    part.buildGraph([0, 2, 6, 9], [2, 1, 2, 4, 3, 0, 3, 1, 0])
elif rank == 1:
    part.buildGraph([0, 5, 8], [2, 5, 1, 7, 4, 1, 3, 7])
elif rank == 2:
    part.buildGraph([0, 3, 6], [3, 7, 6, 5, 7, 8])
elif rank == 3:
    part.buildGraph([0, 5, 7], [4, 3, 5, 6, 8, 7, 6])
part.checkGraph()

meshGraph = code_aster.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh2)
test.assertTrue(meshGraph.debugCheck())
part2 = code_aster.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
outMesh2 = bMesh.applyBalancingStrategy(scotchPart)

checkJoints(outMesh2)

mesh3 = code_aster.IncompleteMesh()
mesh3.readMedFile("petsc04a.mmed")
checkMesh3 = code_aster.Mesh()
checkMesh3.readMedFile("petsc04a.mmed")
test.assertTrue(mesh3.debugCheckFromBaseMesh(checkMesh3))
bMesh = code_aster.MeshBalancer()
bMesh.buildFromBaseMesh(mesh3)
meshGraph = code_aster.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh3)
test.assertTrue(meshGraph.debugCheck())
part2 = code_aster.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
outMesh3 = bMesh.applyBalancingStrategy(scotchPart)

checkJoints(outMesh3)

mesh3 = code_aster.IncompleteMesh()
mesh3.readMedFile("forma02a.mmed")
checkMesh3 = code_aster.Mesh()
checkMesh3.readMedFile("forma02a.mmed")
test.assertEqual(mesh3.debugCheckFromBaseMesh(checkMesh3), True)
bMesh = code_aster.MeshBalancer()
bMesh.buildFromBaseMesh(mesh3)
meshGraph = code_aster.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh3)
part2 = code_aster.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
outMesh3 = bMesh.applyBalancingStrategy(scotchPart)

checkJoints(outMesh3)

FIN()
