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

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

monMaillage = CA.Mesh()
monMaillage.readMedFile("zzzz503a.mmed")

monModel = CA.Model(monMaillage)
monModel.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
monModel.build()


YOUNG = 200000.0
POISSON = 0.3

acier = DEFI_MATERIAU(ELAS=_F(E=YOUNG, NU=POISSON))
# acier.debugPrint(6)
test.assertEqual(acier.getType(), "MATER_SDASTER")

affectMat = CA.MaterialField(monMaillage)
affectMat.addMaterialOnMesh(acier)
affectMat.addMaterialOnGroupOfCells(acier, ["Haut", "Bas"])
affectMat.build()
test.assertEqual(affectMat.getType(), "CHAM_MATER")

imposedDof1 = CA.DisplacementReal()
imposedDof1.setValue(CA.PhysicalQuantityComponent.Dx, 0.0)
imposedDof1.setValue(CA.PhysicalQuantityComponent.Dy, 0.0)
imposedDof1.setValue(CA.PhysicalQuantityComponent.Dz, 0.0)
CharMeca1 = CA.ImposedDisplacementReal(monModel)
CharMeca1.setValue(imposedDof1, "Bas")
CharMeca1.build()
test.assertEqual(CharMeca1.getType(), "CHAR_MECA")

imposedPres1 = CA.PressureReal()
imposedPres1.setValue(CA.PhysicalQuantityComponent.Pres, 1000.0)
CharMeca2 = CA.DistributedPressureReal(monModel)
CharMeca2.setValue(imposedPres1, "Haut")
CharMeca2.build()
test.assertEqual(CharMeca2.getType(), "CHAR_MECA")

# create ant test PhysicalProblem
study = CA.PhysicalProblem(monModel, affectMat)
study.addLoad(CharMeca1)
study.addLoad(CharMeca2)
test.assertTrue(study.computeListOfLoads())
test.assertTrue(study.computeDOFNumbering())
listLoads = study.getListOfLoads()
dofNume = study.getDOFNumbering()
study.computeBehaviourProperty(COMPORTEMENT=(_F(RELATION="VMIS_ISOT_LINE", TOUT="OUI"),))
dofNume = study.getBehaviourProperty()


test.assertEqual(monMaillage.getName(), study.getMesh().getName())
test.assertEqual(monModel.getName(), study.getModel().getName())
test.assertEqual(affectMat.getName(), study.getMaterialField().getName())
test.assertEqual(None, study.getElementaryCharacteristics())

# comparison between DataStructures
test.assertEqual(study.getMesh(), monMaillage)
test.assertEqual(study.getModel(), monModel)
test.assertNotEqual(study.getModel(), monMaillage)
test.assertNotEqual(study.getModel(), "model ?")
test.assertNotEqual(study.getModel(), object())

FIN()
