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

# --------------------------------------------------------------------
# This file is part of
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
# along with   If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

from . import MedFileReader, IncompleteMesh, MeshBalancer, MeshConnectionGraph, PtScotchPartitioner
from . import FieldCharacteristics, SimpleFieldOnNodesReal, Result
from . import SimpleFieldOnCellsReal
from . import MYMED2ASTER_CONNECT, MED_TYPES, ASTER_TYPES


def splitMeshAndFieldsFromMedFile(filename, cellBalancer=False, nodeBalancer=False, outMesh=None):
    """Split a MED mesh and MED fields from a filename

    Arguments:
        filename (str): filename of MED file.
        cellBalancer (bool): True if cell balancer must be return.
        nodeBalancer (bool): True if node balancer must be return.
        outMesh (ParallelMesh): split mesh.

    Returns:
        tuple: first element: split mesh, second element: dict with split fields,
               third element: cell balancer (if asked),
               fourth element: node balancer (if asked).
    """
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


def splitMedFileToResults(filename, fieldToRead, resultType, model=None):
    """Split a MED mesh and MED fields from a filename and return Result

    Arguments:
        filename (str): filename of MED file.
        fieldToRead (dict): dict that matches med field name and aster name (in Result)
        resultType (class): A Result class to instanciate (child of Result)
        model (Model): Model (in case of ELGA field reading)

    Returns:
        Result: results container (type: resultType) with mesh and all readen fields
    """
    if not issubclass(resultType, Result):
        raise TypeError("resultType must be a child class of Result")
    result = resultType()

    # Split mesh and fields
    (mesh, fields) = splitMeshAndFieldsFromMedFile(filename)

    if model is not None:
        mesh = model.getMesh()

    toMedGeoType = dict(zip(ASTER_TYPES, MED_TYPES))

    first = True
    # Loop over field dict given by user
    for medFieldName in fieldToRead:
        # Get aster name of med field
        asterFieldName = fieldToRead[medFieldName]

        # Get characteristics (localization, quantity name) of field from aster name
        fieldChar = FieldCharacteristics(asterFieldName)
        loc = fieldChar.getLocalization()
        qt = fieldChar.getQuantity()
        opt = fieldChar.getOption()
        param = fieldChar.getParameter()

        # Get split field
        curMedFieldDict = fields[medFieldName]
        # Resize output result
        if first:
            result.resize(len(curMedFieldDict))
            first = False
        # Loop over field "time" index
        for index in curMedFieldDict:
            curField = curMedFieldDict[index]
            fieldValues = curField.getValues()
            compName = curField.getComponentName()
            fieldToAdd = None
            # FieldOnNodes case
            if loc == "NOEU":
                sFON = SimpleFieldOnNodesReal(mesh, qt, compName, True)
                nodeNb = mesh.getNumberOfNodes()
                cmpNb = len(compName)
                i = 0
                # Copy values in field
                for iNode in range(nodeNb):
                    for iCmp in range(cmpNb):
                        sFON[iNode, iCmp] = fieldValues[i]
                        i += 1

                # Convert SimpleFieldOnNodes to FieldOnNodes
                fieldToAdd = sFON.toFieldOnNodes()
            # FieldOnCells case
            elif loc in ("ELGA", "ELNO", "ELEM"):
                if model is None:
                    raise NameError("Model is mandatory when ELGA field is asked")
                cumSizes = curField.getCumulatedSizesVector()
                cmpNb = len(compName)
                fieldCmps = curField.getComponentVector()
                maxNbGPNb = -1
                for iCell, curCmpNum in enumerate(fieldCmps):
                    gPNb = curCmpNum / cmpNb
                    if int(gPNb) == gPNb:
                        gPNb = int(gPNb)
                    else:
                        raise NameError("Inconsistent Gauss point number")
                    maxNbGPNb = max(gPNb, maxNbGPNb)
                sFOC = SimpleFieldOnCellsReal(mesh, loc, qt, compName, maxNbGPNb, 1, True)
                # Copy values in field
                if loc == "ELNO":
                    assert len(fieldCmps) == mesh.getNumberOfCells()
                    for iCell, curCmpNum in enumerate(fieldCmps):
                        gPNb = int(curCmpNum / cmpNb)
                        j = 0
                        posInNew = cumSizes[iCell]
                        cellType = mesh.getCellType(iCell)
                        medType = toMedGeoType[cellType]
                        for iPt in range(gPNb):
                            newPos = MYMED2ASTER_CONNECT[medType][iPt]
                            for iCmp in range(cmpNb):
                                sFOC[iCell, iCmp, newPos, 0] = fieldValues[posInNew + j]
                                j += 1
                else:
                    for iCell, curCmpNum in enumerate(fieldCmps):
                        gPNb = int(curCmpNum / cmpNb)
                        j = 0
                        posInNew = cumSizes[iCell]
                        for iPt in range(gPNb):
                            for iCmp in range(cmpNb):
                                sFOC[iCell, iCmp, iPt, 0] = fieldValues[posInNew + j]
                                j += 1
                fED = model.getFiniteElementDescriptor()
                # Convert SimpleFieldOnells to FieldOnCells
                fieldToAdd = sFOC.toFieldOnCells(fED, opt, param)
            else:
                raise NameError("Not yet implemented")
            # Add field to Result
            result.setField(fieldToAdd, asterFieldName, index)

    return result
