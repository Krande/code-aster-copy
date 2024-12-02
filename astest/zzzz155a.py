# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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


from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()


def checkJoints(mesh):
    comm = MPI.ASTER_COMM_WORLD
    l2G = mesh.getLocalToGlobalNodeIds()

    j = 0
    for proc in mesh.getOppositeDomains():
        fJ = mesh.getSendJoint(j)
        gFJ = []
        for i in range(int(len(fJ) / 2)):
            gFJ.append(l2G[fJ[2 * i] - 1])

        sJ = mesh.getReceiveJoint(j)
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


graph = CA.CommGraph()
balancer = CA.ObjectBalancer()
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

bMesh = CA.MeshBalancer()
bMesh2 = CA.MeshBalancer()
# Partitioning from BaseMesh and then from ParallelMesh
# From ParallelMesh, it tests mesh swap in multiple ways
if rank == 0:
    myMesh = CA.Mesh()
    myMesh.readMedFile("fort.20")
    bMesh.buildFromBaseMesh(myMesh)
    outMesh = bMesh.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
    bMesh2.buildFromBaseMesh(outMesh)
    outMesh2 = bMesh2.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
    outMesh3 = bMesh2.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
    outMesh4 = bMesh2.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
elif rank == 1:
    outMesh = bMesh.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
    bMesh2.buildFromBaseMesh(outMesh)
    outMesh2 = bMesh2.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
    outMesh3 = bMesh2.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
    outMesh4 = bMesh2.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
elif rank == 2:
    outMesh = bMesh.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
    bMesh2.buildFromBaseMesh(outMesh)
    outMesh2 = bMesh2.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])
    outMesh3 = bMesh2.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])
    outMesh4 = bMesh2.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
elif rank == 3:
    outMesh = bMesh.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])
    bMesh2.buildFromBaseMesh(outMesh)
    outMesh2 = bMesh2.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
    outMesh3 = bMesh2.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
    outMesh4 = bMesh2.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])

checkMesh = CA.Mesh()
checkMesh.readMedFile("fort.20")

mesh2 = CA.IncompleteMesh()
mesh2.readMedFile("fort.20")
test.assertEqual(mesh2.debugCheckFromBaseMesh(checkMesh), True)

bMesh = CA.MeshBalancer()
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
            [1, 0],
            [0, 2],
            [2, 3],
            [1, 4],
            [4, 5],
            [1, 6],
            [6, 7],
            [0, 8],
            [8, 9],
            [1, 0, 2, 4],
            [4, 2, 3, 5],
            [1, 6, 8, 0],
            [6, 7, 9, 8],
            [1, 4, 10, 6],
            [6, 10, 11, 7],
            [4, 5, 12, 10],
            [10, 12, 13, 11],
            [3, 2, 16, 14],
            [14, 16, 17, 15],
            [2, 0, 8, 16],
            [16, 8, 9, 17],
            [0, 1, 4, 2, 8, 6, 10, 16],
            [8, 6, 10, 16, 9, 7, 11, 17],
            [2, 4, 5, 3, 16, 10, 12, 14],
            [16, 10, 12, 14, 17, 11, 13, 15],
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
            [5, 4],
            [4, 0],
            [0, 1],
            [5, 2],
            [17, 4],
            [2, 3],
            [14, 15],
            [15, 5],
            [16, 17],
            [14, 15, 17, 16],
            [15, 5, 4, 17],
            [14, 6, 7, 15],
            [4, 5, 2, 0],
            [0, 2, 3, 1],
            [10, 12, 13, 11],
            [11, 13, 0, 1],
            [12, 16, 17, 13],
            [13, 17, 4, 0],
            [15, 7, 2, 5],
            [6, 8, 9, 7],
            [7, 9, 3, 2],
            [17, 15, 7, 13, 4, 5, 2, 0],
            [16, 14, 6, 12, 17, 15, 7, 13],
            [12, 6, 8, 10, 13, 7, 9, 11],
            [13, 7, 9, 11, 0, 2, 3, 1],
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
            [16, 17],
            [17, 13],
            [0, 1],
            [1, 13],
            [2, 3],
            [3, 12],
            [14, 15],
            [15, 12],
            [13, 12],
            [4, 6, 7, 5],
            [5, 7, 17, 16],
            [6, 0, 1, 7],
            [7, 1, 13, 17],
            [2, 8, 9, 3],
            [3, 9, 15, 12],
            [8, 10, 11, 9],
            [9, 11, 14, 15],
            [14, 16, 17, 15],
            [15, 17, 13, 12],
            [13, 1, 3, 12],
            [1, 0, 2, 3],
            [10, 4, 6, 8, 11, 5, 7, 9],
            [11, 5, 7, 9, 14, 16, 17, 15],
            [8, 6, 0, 2, 9, 7, 1, 3],
            [9, 7, 1, 3, 15, 17, 13, 12],
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
            [17, 8],
            [8, 9],
            [16, 10],
            [10, 11],
            [15, 17],
            [12, 13],
            [13, 16],
            [17, 16],
            [14, 15],
            [16, 13, 4, 10],
            [10, 4, 5, 11],
            [13, 12, 6, 4],
            [4, 6, 7, 5],
            [14, 15, 2, 0],
            [0, 2, 3, 1],
            [15, 17, 8, 2],
            [2, 8, 9, 3],
            [9, 8, 10, 11],
            [8, 17, 16, 10],
            [14, 12, 13, 15],
            [15, 13, 16, 17],
            [13, 15, 17, 16, 4, 2, 8, 10],
            [4, 2, 8, 10, 5, 3, 9, 11],
            [6, 0, 2, 4, 7, 1, 3, 5],
            [12, 14, 15, 13, 6, 0, 2, 4],
        ],
        True,
    )
checkJoints(outMesh)

part = CA.PtScotchPartitioner()
if rank == 0:
    part.buildGraph([0, 2, 6, 9], [2, 1, 2, 4, 3, 0, 3, 1, 0])
elif rank == 1:
    part.buildGraph([0, 5, 8], [2, 5, 1, 7, 4, 1, 3, 7])
elif rank == 2:
    part.buildGraph([0, 3, 6], [3, 7, 6, 5, 7, 8])
elif rank == 3:
    part.buildGraph([0, 5, 7], [4, 3, 5, 6, 8, 7, 6])
part.checkGraph()

meshGraph = CA.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh2)
test.assertTrue(meshGraph.debugCheck())
part2 = CA.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
outMesh2 = bMesh.applyBalancingStrategy(scotchPart)

checkJoints(outMesh2)

mesh3 = CA.IncompleteMesh()
mesh3.readMedFile("petsc04a.mmed")
checkMesh3 = CA.Mesh()
checkMesh3.readMedFile("petsc04a.mmed")
test.assertTrue(mesh3.debugCheckFromBaseMesh(checkMesh3))
bMesh = CA.MeshBalancer()
bMesh.buildFromBaseMesh(mesh3)
meshGraph = CA.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh3)
test.assertTrue(meshGraph.debugCheck())
part2 = CA.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
outMesh3 = bMesh.applyBalancingStrategy(scotchPart)
innerNodes1 = len(outMesh3.getInnerNodes())

meshGraph = CA.MeshConnectionGraph()
weights = [1] * 55
if rank == 0:
    # Add load on first node
    weights[0] = 10
    meshGraph.buildFromIncompleteMeshWithVertexWeights(mesh3, weights)
elif rank == 1:
    meshGraph.buildFromIncompleteMeshWithVertexWeights(mesh3, weights)
elif rank == 2:
    meshGraph.buildFromIncompleteMeshWithVertexWeights(mesh3, weights)
elif rank == 3:
    meshGraph.buildFromIncompleteMeshWithVertexWeights(mesh3, weights)
test.assertTrue(meshGraph.debugCheck())
part2 = CA.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
outMesh4 = bMesh.applyBalancingStrategy(scotchPart)
index = -1
try:
    index = scotchPart.index(1)
except:
    pass
if index != -1:
    # Node 1 has been loaded
    # The processor that will contain it, will have fewer inner nodes than before
    test.assertTrue(len(outMesh4.getInnerNodes()) < innerNodes1)
else:
    test.assertTrue(len(outMesh4.getInnerNodes()) >= innerNodes1)

checkJoints(outMesh3)

mesh3 = CA.IncompleteMesh()
mesh3.readMedFile("forma02a.mmed")
checkMesh3 = CA.Mesh()
checkMesh3.readMedFile("forma02a.mmed")
test.assertEqual(mesh3.debugCheckFromBaseMesh(checkMesh3), True)
bMesh = CA.MeshBalancer()
bMesh.buildFromBaseMesh(mesh3)
meshGraph = CA.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh3)
part2 = CA.PtScotchPartitioner()
part2.buildGraph(meshGraph)
scotchPart = part2.partitionGraph()
outMesh3 = bMesh.applyBalancingStrategy(scotchPart)

checkJoints(outMesh3)

FIN()
