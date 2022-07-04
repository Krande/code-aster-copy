# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
from code_aster import MPI


code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
pMesh = code_aster.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med"%rank, True)

monModel = code_aster.Model(pMesh)
monModel.addModelingOnMesh(code_aster.Physics.Mechanics,
                              code_aster.Modelings.Tridimensional)
monModel.build()

testMesh = monModel.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

acier = DEFI_MATERIAU(ELAS = _F(E = 2.e11,
                                NU = 0.3,),)

affectMat = code_aster.MaterialField(pMesh)
affectMat.addMaterialOnMesh( acier )
affectMat.build()

testMesh2 = affectMat.getMesh()
test.assertEqual(testMesh2.getType(), "MAILLAGE_P")

charCine = code_aster.MechanicalDirichletBC(monModel)
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dx, 0., "COTE_B")
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dy, 0., "COTE_B")
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dz, 0., "COTE_B")
charCine.build()

charCine2 = code_aster.MechanicalDirichletBC(monModel)
charCine2.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dz, 1., "COTE_H")
charCine2.build()

resu = MECA_STATIQUE(MODELE=monModel, CHAM_MATER=affectMat,
  EXCIT=(_F(CHARGE=charCine), _F(CHARGE=charCine2)),)

resu.printMedFile("fort."+str(rank+40)+".med")

MyFieldOnNodes = resu.getFieldOnNodesReal("DEPL", 1)
sfon = MyFieldOnNodes.exportToSimpleFieldOnNodes()
sfon.build()

val = [0.134202362865, 0.134202362865, 0.154144849556, 0.154144849556]
print(rank, sfon.getValue(4, 1))
test.assertAlmostEqual(sfon.getValue(4, 1), val[rank])

test.printSummary()

FIN()
