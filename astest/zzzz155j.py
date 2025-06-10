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

# This test checks elga field writing in a med file
# with re-reading of this field

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI
import os.path as osp

from code_aster.MedUtils import splitMedFileToResults
from code_aster.Utilities import SharedTmpdir

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

# Mesh reading
mesh1 = CA.ParallelMesh()
mesh1.readMedFile("zzzz155j.med")

# Model definition (modeling melting)
model = AFFE_MODELE(
    MAILLAGE=mesh1,
    AFFE=(
        _F(
            MODELISATION="3D",
            PHENOMENE="MECANIQUE",
            GROUP_MA=(
                "vis",
                "cloison",
                "renfort",
                "VisAppuiTete",
                "CloisonSymetrieX",
                "CloisonSymetrieZ",
                "RenfortSymetrieX",
                "RenfortBloquage",
            ),
        ),
        _F(MODELISATION="DIS_T", PHENOMENE="MECANIQUE", GROUP_MA=("DiscretsContact", "VisBasR")),
    ),
)

VisBasR = mesh1.getCells("VisBasR")

# On cell group VisBasR: 2 Gauss points
gPList = [0] * mesh1.getNumberOfCells()
for iCell in VisBasR:
    gPList[iCell - 1] = 2

# On 3D modeling, 5 or 8 Gauss points
types = mesh1.getMedCellsTypes()
for iCell in range(len(gPList)):
    curType = types[iCell]
    if curType == 310:
        gPList[iCell] = 5
    elif curType == 308:
        gPList[iCell] = 8

compName = ["SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ", "N", "VY", "VZ"]

# Simple field on cells allocation
sFOC = CA.SimpleFieldOnCellsReal(mesh1, "ELGA", "SIEF_R", compName, gPList, 1, True)

# Field values are Gauss point numbers only on current modeling components:
#  - 'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ' for 3D,
#  - 'N', 'VY', 'VZ' for DIS_T
cmpNb = len(compName)
for iCell in range(len(gPList)):
    gPNb = gPList[iCell]
    cmpS = 0
    cmpE = 0
    if gPNb == 2:
        cmpS = 6
        cmpE = 9
    elif gPNb == 5 or gPNb == 8:
        cmpS = 0
        cmpE = 6
    for iPt in range(gPNb):
        for iCmp in range(cmpS, cmpE):
            sFOC.setValue(iCell, iCmp, iPt, 0, gPNb)

fED = model.getFiniteElementDescriptor()
# convert to FieldOnCells
fieldToAdd = sFOC.toFieldOnCells(fED, "RAPH_MECA", "PCONTPR")

result = CA.ElasticResult()
result.resize(1)
result.setField(fieldToAdd, "SIEF_ELGA", 1)
result.setTime(0.0, 1)
result.userName = "MONRESU1"

medfile = ""
with SharedTmpdir("zzzz155j_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu.zzzz155j.med")
    DEFI_FICHIER(UNITE=80, FICHIER=medfile, TYPE="BINARY")
    # Single med file
    IMPR_RESU(
        FICHIER_UNIQUE="OUI", FORMAT="MED", UNITE=80, RESU=_F(RESULTAT=result), VERSION_MED="4.1.0"
    )
    DEFI_FICHIER(ACTION="LIBERER", UNITE=80)

    # Reading of newly created single med file
    mesh2 = CA.ParallelMesh()
    mesh2.readMedFile(medfile)

    model2 = AFFE_MODELE(
        MAILLAGE=mesh2,
        AFFE=(
            _F(
                MODELISATION="3D",
                PHENOMENE="MECANIQUE",
                GROUP_MA=(
                    "vis",
                    "cloison",
                    "renfort",
                    "VisAppuiTete",
                    "CloisonSymetrieX",
                    "CloisonSymetrieZ",
                    "RenfortSymetrieX",
                    "RenfortBloquage",
                ),
            ),
            _F(
                MODELISATION="DIS_T", PHENOMENE="MECANIQUE", GROUP_MA=("DiscretsContact", "VisBasR")
            ),
        ),
    )
    # Split field
    resultCheck = splitMedFileToResults(
        medfile, {"MONRESU1SIEF_ELGA": "SIEF_ELGA"}, CA.ElasticResult, model2
    )

    fieldCheck = resultCheck.getField("SIEF_ELGA", 1)
    sfieldCheck = fieldCheck.toSimpleFieldOnCells()

    VisBasR = mesh2.getCells("VisBasR")
    setVisBarR = set(VisBasR)

    nbCells = sfieldCheck.getNumberOfCells()
    nbCmp = sfieldCheck.getNumberOfComponents()
    for iCell in range(nbCells):
        nbPt = sfieldCheck.getNumberOfPointsOfCell(iCell)
        cmpS = 0
        cmpE = 0
        # Only on VisBasR: 2 Gauss points
        if nbPt == 2 and iCell + 1 in setVisBarR:
            cmpS = 6
            cmpE = 9
        elif nbPt == 5 or nbPt == 8:
            cmpS = 0
            cmpE = 6
        for iCmp in range(cmpS, cmpE):
            for iPt in range(nbPt):
                test.assertEqual(nbPt, sfieldCheck.getValue(iCell, iCmp, iPt, 0))

FIN()
