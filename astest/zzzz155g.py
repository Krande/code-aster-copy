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

from code_aster.Commands import *
from code_aster import CA
from code_aster.MedUtils import splitMeshAndFieldsFromMedFile


# CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))
DEBUT(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

filename = "zzzz155d.med"

(MecaMail2021, fields) = splitMeshAndFieldsFromMedFile(filename)
# fieldToRead = {
# "TEMP____TEMP": "TEMP",
# "IRRA____IRRA": "IRRA",
# }
fieldToRead = {
    "00000008DEPL": "DEPL"
}  # "00000008SIEF_ELGA": "SIEF_ELGA", "00000008SIEF_ELNO": "SIEF_ELNO"}
tmp = None
# results = [TEMPER, IRRAD]
results = CA.NonLinearResult()
# Loop over field dict given by user
for idx, medFieldName in enumerate(fieldToRead):
    # Get aster name of med field
    asterFieldName = fieldToRead[medFieldName]

    # Get characteristics (localization, quantity name) of field from aster name
    fieldChar = CA.FieldCharacteristics(asterFieldName)
    loc = fieldChar.getLocalization()
    qt = fieldChar.getQuantity()
    opt = fieldChar.getOption()
    param = fieldChar.getParameter()
    del fieldChar

    # Get split field
    curMedFieldDict = fields[medFieldName]
    # Resize output result
    # results[idx].resize(len(curMedFieldDict))
    results.resize(len(curMedFieldDict))
    # Loop over field "time" index
    for index, curTime in curMedFieldDict["id2time"]:
        curField = curMedFieldDict[index]
        tmp = curField
        # print("LA", type(tmp))
        fieldValues = curField.getValues()
        compName = curField.getComponentName()
        fieldToAdd = None
        # FieldOnNodes case
        if loc == "NOEU":
            sFON = CA.SimpleFieldOnNodesReal(MecaMail2021, qt, compName, True)
            nodeNb = MecaMail2021.getNumberOfNodes()
            cmpNb = len(compName)
            i = 0
            # Copy values in field
            for iNode in range(nodeNb):
                for iCmp in range(cmpNb):
                    sFON[iNode, iCmp] = fieldValues[i]
                    i += 1

            # Convert SimpleFieldOnNodes to FieldOnNodes
            fieldToAdd = sFON.toFieldOnNodes()
        # Add field to Result
        # results[idx].setField(fieldToAdd, asterFieldName, index)
        # results[idx].setTime(curTime, index)
        results.setField(fieldToAdd, asterFieldName, index)
        results.setTime(curTime, index)

values = curField.getValues()
print(values)
test.assertEqual(len(values), 81)
test.assertAlmostEqual(values[-3], 0.0)
test.assertAlmostEqual(values[-2], -0.1153846)
test.assertAlmostEqual(values[-1], 0.5)

FIN(PROC0="NON")
