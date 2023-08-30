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

from ...Objects import MedFileReader, IncompleteMesh, MeshBalancer, MeshConnectionGraph
from ...Objects import PtScotchPartitioner


def splitMeshAndFieldsFromMedFile(filename, cellBalancer=False, nodeBalancer=False, outMesh=None):
    fr = MedFileReader()
    fr.openParallel(filename)
    mesh = IncompleteMesh()
    mesh.readMedFile(filename)
    bMesh = MeshBalancer()
    bMesh.buildFromBaseMesh(mesh)
    meshGraph = MeshConnectionGraph()
    meshGraph.buildFromIncompleteMesh(mesh)

    part2 = PtScotchPartitioner()
    part2.buildGraph(meshGraph)
    scotchPart = part2.partitionGraph()

    outMesh = bMesh.applyBalancingStrategy(scotchPart, outMesh)

    medCellTypes = outMesh.getMedCellsTypes()
    oldType = -999
    sortedTypes = []
    checkDict = {}
    for type in medCellTypes:
        if oldType != type:
            sortedTypes.append(type)
            oldType = type
            if checkDict.get(type) is not None:
                raise RuntimeError("Cells are not sorted by type")
            checkDict[type] = 1

    nBalancer = bMesh.getNodeObjectBalancer()
    cBalancer = bMesh.getCellObjectBalancer()

    nbField = fr.getFieldNumber()
    fieldDict = {}

    for i in range(nbField):
        curField = fr.getField(i)
        nbSeq = curField.getSequenceNumber()
        fieldDict[curField.getName()] = {}
        for j in range(nbSeq):
            curSeq = curField.getSequence(j)
            assert curSeq[0] == curSeq[1]
            supportEnt = curField.getAllSupportEntitiesAtSequence(curSeq[0], curSeq[1])
            if supportEnt[0] == 3:
                valuesVec = curField.getValuesAtSequenceOnNodes(curSeq[0], curSeq[1])
                out = nBalancer.balanceMedVectorOverProcessesWithRenumbering(valuesVec)
                fieldDict[curField.getName()][curSeq[0]] = out
            else:
                valuesVec = curField.getValuesAtSequenceOnCellTypesList(
                    curSeq[0], curSeq[1], sortedTypes
                )
                out = cBalancer.balanceMedVectorOverProcessesWithRenumbering(valuesVec)
                fieldDict[curField.getName()][curSeq[0]] = out
    toReturn = (outMesh, fieldDict)
    if cellBalancer:
        toReturn += (cBalancer,)
    if nodeBalancer:
        toReturn += (nBalancer,)
    return toReturn
