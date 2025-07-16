# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
from code_aster.MedUtils import (
    convertMedFieldToAster,
    splitMedFileToResults,
    splitMeshAndFieldsFromMedFile,
)

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()


def splitEntitySet(nbElemT, rank, nbProcs):
    nbElemL = int(nbElemT / nbProcs)
    start = rank * nbElemL + 1
    if rank == nbProcs - 1:
        end = nbProcs * nbElemL
        nbElemL = nbElemL - (end - nbElemT)
    return (nbElemL, start)


filename = "zzzz155d.med"

myTuple = splitMeshAndFieldsFromMedFile(filename, True)

# Small example of what can be done to change component names while reading med file
fieldList = myTuple[1]
curMedDepl = fieldList["00000008DEPL"][1]
# Displacement field is now a temperature field
fieldToAdd = convertMedFieldToAster(
    curMedDepl, "TEMP", myTuple[0], ["TEMP", "TEMP_MIL", "TEMP_INF"]
)
# add in a ThermalResult
thermRes = CA.ThermalResult()
thermRes.resize(1)
thermRes.setField(fieldToAdd, "TEMP", 1)
thermRes.setTime(1.0, 1)

valuesDeplRef = [
    0.000000,
    0.000000,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    0.173077,
    0.115385,
    0.500000,
    0.000000,
    0.000000,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    -0.173077,
    0.115385,
    0.500000,
    0.000000,
    0.000000,
    1.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    1.000000,
    0.173077,
    -0.115385,
    0.500000,
    -0.173077,
    -0.115385,
    0.500000,
    0.000000,
    0.115385,
    0.500000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    0.000000,
    1.000000,
    0.000000,
    -0.115385,
    0.500000,
]

profileDepl = {
    0: 0,
    1: 1,
    4: 2,
    5: 3,
    8: 4,
    9: 5,
    11: 6,
    12: 7,
    13: 8,
    15: 9,
    16: 10,
    17: 11,
    20: 12,
    21: 13,
    22: 14,
    24: 15,
    25: 16,
    26: 17,
}

splitMesh = myTuple[0]
loc2Glob = splitMesh.getLocalToGlobalNodeIds()
depl = myTuple[1]["00000008DEPL"][1]
valuesDepl = depl.getValues()
for count, i in enumerate(loc2Glob):
    if i in profileDepl:
        posInRef = profileDepl[i]
        # 3 components: DX, DY, DZ
        for j in range(3):
            test.assertAlmostEqual(
                abs(valuesDepl[3 * count + j] - valuesDeplRef[posInRef * 3 + j]), 0, delta=1e-6
            )

valuesSiefRef = [
    88657958936.584000,
    84719063516.309814,
    262013106735.868134,
    -0.000008,
    -22050236751.245529,
    -10761262413.889463,
    88657958936.584030,
    84719063516.309845,
    262013106735.868195,
    -0.000004,
    -5908343130.411303,
    -10761262413.889467,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    -0.000006,
    -22050236751.245514,
    10761262413.889471,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    -0.000003,
    -5908343130.411303,
    10761262413.889450,
    -122455978.004181,
    -14822613812.167831,
    205516479062.948395,
    -0.000027,
    -22050236751.245548,
    -10761262413.889454,
    -122455978.004135,
    -14822613812.167786,
    205516479062.948456,
    -0.000012,
    -5908343130.411332,
    -10761262413.889458,
    -122455978.004120,
    -14822613812.167801,
    205516479062.948425,
    -0.000027,
    -22050236751.245522,
    10761262413.889481,
    -122455978.004135,
    -14822613812.167816,
    205516479062.948425,
    -0.000016,
    -5908343130.411311,
    10761262413.889484,
    88657958936.584045,
    84719063516.309845,
    262013106735.868195,
    -0.000004,
    5908343130.411306,
    -10761262413.889450,
    88657958936.584000,
    84719063516.309814,
    262013106735.868073,
    -0.000004,
    22050236751.245514,
    -10761262413.889431,
    88657958936.584045,
    84719063516.309845,
    262013106735.868195,
    -0.000001,
    5908343130.411303,
    10761262413.889494,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    0.000000,
    22050236751.245502,
    10761262413.889503,
    -122455978.004105,
    -14822613812.167786,
    205516479062.948456,
    -0.000015,
    5908343130.411319,
    -10761262413.889471,
    -122455978.004135,
    -14822613812.167816,
    205516479062.948364,
    -0.000010,
    22050236751.245525,
    -10761262413.889467,
    -122455978.004135,
    -14822613812.167801,
    205516479062.948456,
    -0.000009,
    5908343130.411338,
    10761262413.889481,
    -122455978.004150,
    -14822613812.167801,
    205516479062.948395,
    -0.000002,
    22050236751.245548,
    10761262413.889481,
    -122455978.004135,
    -14822613812.167786,
    205516479062.948456,
    -0.000027,
    22050236751.245495,
    10761262413.889458,
    -122455978.004166,
    -14822613812.167847,
    205516479062.948395,
    -0.000011,
    5908343130.411278,
    10761262413.889458,
    -122455978.004089,
    -14822613812.167770,
    205516479062.948486,
    -0.000025,
    22050236751.245476,
    -10761262413.889488,
    -122455978.004120,
    -14822613812.167801,
    205516479062.948425,
    -0.000011,
    5908343130.411266,
    -10761262413.889494,
    88657958936.584045,
    84719063516.309845,
    262013106735.868195,
    -0.000007,
    22050236751.245518,
    10761262413.889458,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    -0.000003,
    5908343130.411300,
    10761262413.889463,
    88657958936.584045,
    84719063516.309845,
    262013106735.868195,
    -0.000008,
    22050236751.245502,
    -10761262413.889492,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    -0.000004,
    5908343130.411291,
    -10761262413.889488,
    -122455978.004166,
    -14822613812.167831,
    205516479062.948395,
    -0.000013,
    -5908343130.411266,
    10761262413.889458,
    -122455978.004120,
    -14822613812.167770,
    205516479062.948456,
    -0.000011,
    -22050236751.245480,
    10761262413.889454,
    -122455978.004166,
    -14822613812.167831,
    205516479062.948395,
    -0.000004,
    -5908343130.411293,
    -10761262413.889490,
    -122455978.004135,
    -14822613812.167801,
    205516479062.948425,
    0.000000,
    -22050236751.245502,
    -10761262413.889494,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    -0.000003,
    -5908343130.411281,
    10761262413.889463,
    88657958936.584045,
    84719063516.309845,
    262013106735.868195,
    -0.000002,
    -22050236751.245499,
    10761262413.889458,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    -0.000002,
    -5908343130.411293,
    -10761262413.889490,
    88657958936.584015,
    84719063516.309814,
    262013106735.868134,
    -0.000001,
    -22050236751.245502,
    -10761262413.889492,
]

# +24 because there is 24 seg2 before hexa8
# Element global numbering is made of: 1, ..., 24 seg2 et 25, ..., 32 hexa8
profileSief = {1 + 24: 0, 2 + 24: 1, 5 + 24: 2, 6 + 24: 3}

# Split cells like if they were read by MedFileReader
split = splitEntitySet(32, rank, size)
if rank == 0:
    cellGlobalId = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 25, 26, 27, 28]
elif rank == 1:
    cellGlobalId = [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 29, 30, 31, 32]
else:
    raise NameError("Test only on 2 procs")

cBalancer = myTuple[2]
# Balance cell global id over processes in order to be able to compare ELGA field
bCellGlobId = cBalancer.balanceVectorOverProcesses(cellGlobalId)
renumber = cBalancer.getRenumbering()
bCellGlobIdR = [0] * len(bCellGlobId)
for i in range(len(bCellGlobId)):
    bCellGlobIdR[renumber[i] - 1] = bCellGlobId[i]

sief = myTuple[1]["00000008SIEF_ELGA"][1]
cumSizes = sief.getCumulatedSizesVector()
valuesSief = sief.getValues()

nbCmp = 6 * 8
for i in range(len(bCellGlobId)):
    numGlob = bCellGlobIdR[i]
    if numGlob in profileSief:
        posInRef = profileSief[numGlob]
        posInNew = cumSizes[i]
        for j in range(nbCmp):
            test.assertAlmostEqual(
                abs(valuesSief[posInNew + j] - valuesSiefRef[(posInRef) * nbCmp + j]), 0, delta=1e-6
            )

deplChar = CA.FieldCharacteristics("DEPL")
deplQt = deplChar.getQuantity()
sFON = CA.SimpleFieldOnNodesReal(splitMesh, deplQt, depl.getComponentName(), True)

nodeNb = splitMesh.getNumberOfNodes()
cmpNb = len(depl.getComponentName())
i = 0
for iNode in range(nodeNb):
    for iCmp in range(cmpNb):
        sFON[iNode, iCmp] = valuesDepl[i]
        i += 1

fON = sFON.toFieldOnNodes()

model = AFFE_MODELE(
    MAILLAGE=splitMesh, AFFE=_F(GROUP_MA="DEMICUBE", PHENOMENE="MECANIQUE", MODELISATION="3D")
)

siefChar = CA.FieldCharacteristics("SIEF_ELGA")
siefQt = siefChar.getQuantity()
siefLoc = siefChar.getLocalization()
sFOC = CA.SimpleFieldOnCellsReal(splitMesh, siefLoc, siefQt, sief.getComponentName(), 8, 1, True)

cmpNb = len(sief.getComponentName())
siefCmps = sief.getComponentVector()
for iCell, curCmpNum in enumerate(siefCmps):
    posInSief = cumSizes[iCell]
    gPNb = curCmpNum / cmpNb
    if int(gPNb) == gPNb:
        gPNb = int(gPNb)
    else:
        raise NameError("Inconsistent Gauss point number")
    j = 0
    posInNew = cumSizes[iCell]
    for iPt in range(gPNb):
        for iCmp in range(cmpNb):
            sFOC[iCell, iCmp, iPt, 1] = valuesSief[posInNew + j]
            j += 1
fED = model.getFiniteElementDescriptor()
fOC = sFOC.toFieldOnCells(fED, "TOU_INI_ELGA", " ")

result = CA.NonLinearResult()
result.resize(1)
result.setField(fON, "DEPL", 1)
result.setField(fOC, "SIEF_ELGA", 1)

myResu = splitMedFileToResults(
    filename,
    {"00000008DEPL": "DEPL", "00000008SIEF_ELGA": "SIEF_ELGA", "00000008SIEF_ELNO": "SIEF_ELNO"},
    CA.NonLinearResult,
    model,
)

deplSplit = myResu.getField("DEPL", 1)
siefSplit = myResu.getField("SIEF_ELGA", 1)
siefSplit2 = myResu.getField("SIEF_ELNO", 1)

MA = CA.ParallelMesh()
MA.readMedFile("zzzz155d/" + str(rank) + ".med", partitioned=True)
MO = AFFE_MODELE(
    MAILLAGE=MA, AFFE=_F(GROUP_MA="DEMICUBE", PHENOMENE="MECANIQUE", MODELISATION="3D")
)
MAT = DEFI_MATERIAU(ELAS=_F(E=2.1e11, NU=0.3))
CHAM_MAT = AFFE_MATERIAU(MAILLAGE=MA, AFFE=_F(TOUT="OUI", MATER=MAT))

DEFI_FICHIER(UNITE=81, FICHIER="zzzz155d/" + str(rank) + ".med", TYPE="BINARY", ACCES="OLD")
MEC = LIRE_RESU(
    FORMAT="MED",
    MODELE=MO,
    TYPE_RESU="EVOL_ELAS",
    TOUT_ORDRE="OUI",
    UNITE=81,
    FORMAT_MED=(
        _F(NOM_CHAM_MED="0000000aDEPL", NOM_CHAM="DEPL"),
        _F(NOM_CHAM_MED="0000000aSIEF_ELGA", NOM_CHAM="SIEF_ELGA"),
        _F(NOM_CHAM_MED="0000000aSIEF_ELNO", NOM_CHAM="SIEF_ELNO"),
    ),
)
DEFI_FICHIER(ACTION="LIBERER", UNITE=81)

deplSplitRef = MEC.getField("DEPL", 1)
siefSplitRef = MEC.getField("SIEF_ELGA", 1)
siefSplitRef2 = MEC.getField("SIEF_ELNO", 1)

import numpy as np

tab1 = np.array(siefSplit.getValues())
tab2 = np.array(siefSplitRef.getValues())

# La comparaison se fait à 10^-4 près car les valeurs de sief_elga sont de l'ordre de 10^10
if np.allclose(tab1, tab2, atol=1e-4) is not True:
    test.assertTrue(False)
else:
    test.assertTrue(True)

tab1 = np.array(siefSplit2.getValues())
tab2 = np.array(siefSplitRef2.getValues())
if np.allclose(tab1, tab2, atol=1e-4) is not True:
    test.assertTrue(False)
else:
    test.assertTrue(True)

test.assertAlmostEqual(myResu.getTime(1), MEC.getTime(1))

FIN()
