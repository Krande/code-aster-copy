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


mesh2 = CA.IncompleteMesh()
mesh2.readMedFile("fort.20")

meshGraph = CA.MeshConnectionGraph()
meshGraph.buildFromIncompleteMesh(mesh2)
test.assertTrue(meshGraph.debugCheck())

part2 = CA.ParMetisPartitioner()
part2.buildGraph(meshGraph)
metisPart = part2.partitionGraph()

bMesh = CA.MeshBalancer()
bMesh.buildFromBaseMesh(mesh2)
outMesh2 = bMesh.applyBalancingStrategy(metisPart)

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
part2 = CA.ParMetisPartitioner()
part2.buildGraph(meshGraph)
metisPart = part2.partitionGraph()
outMesh3 = bMesh.applyBalancingStrategy(metisPart)

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
part2 = CA.ParMetisPartitioner()
part2.buildGraph(meshGraph)
metisPart = part2.partitionGraph()
outMesh3 = bMesh.applyBalancingStrategy(metisPart)

checkJoints(outMesh3)

FIN()
