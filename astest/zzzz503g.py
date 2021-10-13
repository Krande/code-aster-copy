# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

import code_aster
from code_aster.Commands import *
import libaster
import numpy as np

code_aster.init("--test")

test = code_aster.TestCase()

monMaillage = code_aster.Mesh()
monMaillage.readMedFile("zzzz503a.mmed")

monModel = code_aster.Model(monMaillage)
monModel.addModelingOnMesh(code_aster.Physics.Mechanics, code_aster.Modelings.Tridimensional)
monModel.build()


YOUNG = 200000.0
POISSON = 0.3

acier = DEFI_MATERIAU(ELAS=_F(E=YOUNG,
                              NU=POISSON,),)
# acier.debugPrint(6)
test.assertEqual(acier.getType(), "MATER_SDASTER")

affectMat = code_aster.MaterialField(monMaillage)
affectMat.addMaterialsOnMesh(acier)
affectMat.addMaterialsOnGroupOfCells(acier, ['Haut', 'Bas'])
affectMat.buildWithoutExternalStateVariables()
test.assertEqual(affectMat.getType(), "CHAM_MATER")

imposedDof1 = code_aster.DisplacementReal()
imposedDof1.setValue(code_aster.PhysicalQuantityComponent.Dx, 0.0)
imposedDof1.setValue(code_aster.PhysicalQuantityComponent.Dy, 0.0)
imposedDof1.setValue(code_aster.PhysicalQuantityComponent.Dz, 0.0)
CharMeca1 = code_aster.ImposedDisplacementReal(monModel)
CharMeca1.setValue(imposedDof1, "Bas")
CharMeca1.build()
test.assertEqual(CharMeca1.getType(), "CHAR_MECA")

imposedPres1 = code_aster.PressureReal()
imposedPres1.setValue(code_aster.PhysicalQuantityComponent.Pres, 1000.)
CharMeca2 = code_aster.DistributedPressureReal(monModel)
CharMeca2.setValue(imposedPres1, "Haut")
CharMeca2.build()
test.assertEqual(CharMeca2.getType(), "CHAR_MECA")

# create ant test PhysicalProblem
study = code_aster.PhysicalProblem(monModel, affectMat)
study.addLoad(CharMeca1)
study.addLoad(CharMeca2)
test.assertTrue(study.computeListOfLoads())
test.assertTrue(study.computeDOFNumbering())
listLoads = study.getListOfLoads()
dofNume = study.getDOFNumbering()
dofNume = study.getBehaviourProperty()

test.assertEqual(monMaillage.getName(), study.getMesh().getName())
test.assertEqual(monModel.getName(), study.getModel().getName())
test.assertEqual(affectMat.getName(), study.getMaterialField().getName())
test.assertEqual(None, study.getElementaryCharacteristics())


FIN()
