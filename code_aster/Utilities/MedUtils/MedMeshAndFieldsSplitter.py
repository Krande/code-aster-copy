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

from enum import Enum, auto

from . import MedFileReader, IncompleteMesh, MeshBalancer, MeshConnectionGraph, PtScotchPartitioner
from . import MedFileAccessType, ParallelMesh
from . import FieldCharacteristics, SimpleFieldOnNodesReal, Result
from . import SimpleFieldOnCellsReal
from . import MYMED2ASTER_CONNECT, MED_TYPES, ASTER_TYPES
from . import toAsterGeoType


class MedSequenceToOrder(Enum):
    Both = auto()
    TimeStep = auto()
    Iteration = auto()


def _splitMeshAndFieldsFromMedFile(
    filename,
    cellBalancer=False,
    nodeBalancer=False,
    outMesh=None,
    nodeGrpToGather=[],
    deterministic=False,
    parallel=True,
    sequenceToOrder=MedSequenceToOrder.Both,
):
    """Split a MED mesh and MED fields from a filename

    Arguments:
        filename (Path|str): filename of MED file.
        cellBalancer (bool): True if cell balancer must be return.
        nodeBalancer (bool): True if node balancer must be return.
        outMesh (ParallelMesh): split mesh.
        nodeGrpToGather (list of list): list of list of node groups to gather
            each item of list will be gather on a sigle MPI processor
        parallel (bool): True if med file must be open in parallel
        sequenceToOrder (MedSequenceToOrder): strategy to convert med sequence to aster order od

    Returns:
        tuple: first element: split mesh, second element: dict with split fields,
               third element: cell balancer (if asked),
               fourth element: node balancer (if asked).
    """
    fr = MedFileReader()
    nBalancer = None
    cBalancer = None

    if parallel:
        fr.openParallel(filename, MedFileAccessType.MedReadOnly)
        mesh = IncompleteMesh()
        mesh.readMedFile(filename, verbose=0)
        idsToGather = []
        if nodeGrpToGather != []:
            for tab1 in nodeGrpToGather:
                curIdsToGather = []
                for grpName in tab1:
                    curIdsToGather += mesh.getNodesFromGroup(grpName)
                idsToGather.append(curIdsToGather)
        bMesh = MeshBalancer()
        bMesh.buildFromBaseMesh(mesh)
        meshGraph = MeshConnectionGraph()
        meshGraph.buildFromIncompleteMesh(mesh)

        part2 = PtScotchPartitioner()
        part2.buildGraph(meshGraph, idsToGather)
        scotchPart = part2.partitionGraph(deterministic)

        outMesh = bMesh.applyBalancingStrategy(scotchPart, outMesh)

        nBalancer = bMesh.getNodeObjectBalancer()
        cBalancer = bMesh.getCellObjectBalancer()
    else:
        fr.open(filename, MedFileAccessType.MedReadOnly)
        outMesh = ParallelMesh()
        outMesh.readMedFile(filename, partitioned=True, verbose=0)

    cellTypes = outMesh.getAllMedCellsTypes()
    asterTypes = [toAsterGeoType(item) for item in cellTypes]
    cellTypes = [x for _, x in sorted(zip(asterTypes, cellTypes))]

    nbField = fr.getFieldNumber()
    fieldDict = {}

    for i in range(nbField):
        curField = fr.getField(i)
        nbSeq = curField.getSequenceNumber()
        fieldDict[curField.getName()] = {}
        fieldDict[curField.getName()]["id2time"] = []
        for j in range(nbSeq):
            curSeq = curField.getSequence(j)
            if sequenceToOrder == MedSequenceToOrder.Both:
                if curSeq[0] != curSeq[1]:
                    raise NameError(
                        "Time step and iteration number must be equal in med field\n"
                        "Try to change the sequenceToOrder argument"
                    )
                orderNb = curSeq[1]
            elif sequenceToOrder == MedSequenceToOrder.TimeStep:
                orderNb = curSeq[0]
            elif sequenceToOrder == MedSequenceToOrder.Iteration:
                orderNb = curSeq[1]
            curTime = curField.getTime(j)
            fieldDict[curField.getName()]["id2time"].append((orderNb, curTime))
            supportEnt = curField.getAllSupportEntitiesAtSequence(curSeq[0], curSeq[1])
            if supportEnt[0] == 3:
                valuesVec = curField.getValuesAtSequenceOnNodes(curSeq[0], curSeq[1])
                if parallel:
                    out = nBalancer.balanceMedVectorOverProcessesWithRenumbering(valuesVec)
                else:
                    out = valuesVec
                fieldDict[curField.getName()][orderNb] = out
            else:
                valuesVec = curField.getValuesAtSequenceOnCellTypesList(
                    curSeq[0], curSeq[1], cellTypes
                )
                if parallel:
                    out = cBalancer.balanceMedVectorOverProcessesWithRenumbering(valuesVec)
                else:
                    out = valuesVec
                fieldDict[curField.getName()][orderNb] = out
    toReturn = (outMesh, fieldDict)
    if cellBalancer:
        toReturn += (cBalancer,)
    if nodeBalancer:
        toReturn += (nBalancer,)
    return toReturn


def splitMeshAndFieldsFromMedFile(
    filename,
    cellBalancer=False,
    nodeBalancer=False,
    outMesh=None,
    nodeGrpToGather=[],
    deterministic=False,
    sequenceToOrder=MedSequenceToOrder.Both,
):
    """Split a MED mesh and MED fields from a filename

    Arguments:
        filename (Path|str): filename of MED file.
        cellBalancer (bool): True if cell balancer must be return.
        nodeBalancer (bool): True if node balancer must be return.
        outMesh (ParallelMesh): split mesh.
        nodeGrpToGather (list of list): list of list of node groups to gather
            each item of list will be gather on a sigle MPI processor
        sequenceToOrder (MedSequenceToOrder): strategy to convert med sequence to aster order od

    Returns:
        tuple: first element: split mesh, second element: dict with split fields,
               third element: cell balancer (if asked),
               fourth element: node balancer (if asked).
    """
    return _splitMeshAndFieldsFromMedFile(
        filename,
        cellBalancer,
        nodeBalancer,
        outMesh,
        nodeGrpToGather,
        deterministic,
        parallel=True,
        sequenceToOrder=sequenceToOrder,
    )


def _splitMedFileToResults(
    filename,
    fieldToRead,
    resultType,
    model=None,
    parallel=True,
    sequenceToOrder=MedSequenceToOrder.Both,
):
    """Split a MED mesh and MED fields from a filename and return Result

    Arguments:
        filename (Path|str): filename of MED file.
        fieldToRead (dict): dict that matches med field name and aster name (in Result)
        resultType (class): A Result class to instanciate (child of Result)
        model (Model): Model (in case of ELGA field reading)
        parallel (bool): True if med file must be open in parallel
        sequenceToOrder (MedSequenceToOrder): strategy to convert med sequence to aster order od

    Returns:
        Result: results container (type: resultType) with mesh and all readen fields
    """
    if not issubclass(resultType, Result):
        raise TypeError("resultType must be a child class of Result")
    result = resultType()

    # Split mesh and fields
    (mesh, fields) = _splitMeshAndFieldsFromMedFile(
        filename, parallel=parallel, sequenceToOrder=sequenceToOrder
    )

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
        if fields.get(medFieldName) is None:
            fieldListStr = ", ".join(list(fields.keys()))
            raise NameError(
                "Field " + medFieldName + " is missing. Fields in file: " + fieldListStr
            )
        curMedFieldDict = fields[medFieldName]
        # Resize output result
        if first:
            result.resize(len(curMedFieldDict))
            first = False
        # Loop over field "time" index
        for index, curTime in curMedFieldDict["id2time"]:
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
                gPList = []
                for iCell, curCmpNum in enumerate(fieldCmps):
                    gPNb = curCmpNum / cmpNb
                    if int(gPNb) == gPNb:
                        gPNb = int(gPNb)
                    else:
                        raise NameError("Inconsistent Gauss point number")
                    gPList.append(gPNb)
                sFOC = SimpleFieldOnCellsReal(mesh, loc, qt, compName, gPList, 1, True)
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
                                sFOC.setValue(iCell, iCmp, newPos, 0, fieldValues[posInNew + j])
                                j += 1
                else:
                    for iCell, curCmpNum in enumerate(fieldCmps):
                        gPNb = int(curCmpNum / cmpNb)
                        j = 0
                        posInNew = cumSizes[iCell]
                        for iPt in range(gPNb):
                            for iCmp in range(cmpNb):
                                sFOC.setValue(iCell, iCmp, iPt, 0, fieldValues[posInNew + j])
                                j += 1
                fED = model.getFiniteElementDescriptor()
                # Convert SimpleFieldOnells to FieldOnCells
                fieldToAdd = sFOC.toFieldOnCells(fED, opt, param)
            else:
                raise NameError("Not yet implemented")
            # Add field to Result
            result.setField(fieldToAdd, asterFieldName, index)
            result.setTime(curTime, index)

    return result


def splitMedFileToResults(
    filename, fieldToRead, resultType, model=None, sequenceToOrder=MedSequenceToOrder.Both
):
    """Split a MED mesh and MED fields from a filename and return Result

    Arguments:
        filename (Path|str): filename of MED file.
        fieldToRead (dict): dict that matches med field name and aster name (in Result)
        resultType (class): A Result class to instanciate (child of Result)
        model (Model): Model (in case of ELGA field reading)
        sequenceToOrder (MedSequenceToOrder): strategy to convert med sequence to aster order od

    Returns:
        Result: results container (type: resultType) with mesh and all readen fields
    """
    return _splitMedFileToResults(
        filename, fieldToRead, resultType, model, parallel=True, sequenceToOrder=sequenceToOrder
    )


def readMedFileToResults(
    filename, fieldToRead, resultType, model=None, sequenceToOrder=MedSequenceToOrder.Both
):
    """Read a MED mesh and MED fields from a filename and return Result

    Arguments:
        filename (Path|str): filename of MED file.
        fieldToRead (dict): dict that matches med field name and aster name (in Result)
        resultType (class): A Result class to instanciate (child of Result)
        model (Model): Model (in case of ELGA field reading)
        sequenceToOrder (MedSequenceToOrder): strategy to convert med sequence to aster order od

    Returns:
        Result: results container (type: resultType) with mesh and all readen fields
    """
    return _splitMedFileToResults(
        filename, fieldToRead, resultType, model, parallel=False, sequenceToOrder=sequenceToOrder
    )
