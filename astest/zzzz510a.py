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

import numpy as N
import code_aster
from code_aster.Commands import *
from code_aster import MPI

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

graph = code_aster.CommGraph()
balancer = code_aster.ObjectBalancer()
a = [i + rank * 10 for i in range(10)]
print(rank, a)
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

import os

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

mesh2 = code_aster.IncompleteMesh()
mesh2.readMedFile("fort.20")

bMesh = code_aster.MeshBalancer()
bMesh.buildFromBaseMesh(mesh2)
if rank == 0:
    outMesh = bMesh.applyBalancingStrategy([1, 2, 9, 11, 17, 19, 25, 31])
elif rank == 1:
    outMesh = bMesh.applyBalancingStrategy([5, 6, 13, 15, 18, 20, 26, 32])
elif rank == 2:
    outMesh = bMesh.applyBalancingStrategy([7, 8, 14, 16, 22, 24, 28, 30])
elif rank == 3:
    outMesh = bMesh.applyBalancingStrategy([3, 4, 10, 12, 21, 23, 27, 29])

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
part.writeGraph("/home/H85256/dev/codeaster/tmp/Bis" + str(rank) + ".graph")
# print(part.partitionGraph())

FIN()
