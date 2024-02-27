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


from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI


CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
pMesh = CA.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med" % rank, partitioned=True)

monModel = CA.Model(pMesh)
monModel.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
monModel.build()

testMesh = monModel.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

acier = DEFI_MATERIAU(ELAS=_F(E=2.0e11, NU=0.3))

affectMat = CA.MaterialField(pMesh)
affectMat.addMaterialOnMesh(acier)
affectMat.build()

testMesh2 = affectMat.getMesh()
test.assertEqual(testMesh2.getType(), "MAILLAGE_P")

charCine = CA.MechanicalDirichletBC(monModel)
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dx, 0.0, "COTE_B")
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dy, 0.0, "COTE_B")
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dz, 0.0, "COTE_B")
charCine.build()

charMeca = AFFE_CHAR_MECA(
    MODELE=monModel, DOUBLE_LAGRANGE="NON", DDL_IMPO=_F(GROUP_NO=("COTE_H"), DZ=1.0)
)

resu = MECA_STATIQUE(
    MODELE=monModel, CHAM_MATER=affectMat, EXCIT=(_F(CHARGE=charCine), _F(CHARGE=charMeca))
)

resu.printMedFile("fort." + str(rank + 40) + ".med")

MyFieldOnNodes = resu.getField("DEPL", 1)
sfon = MyFieldOnNodes.toSimpleFieldOnNodes()
sfon.build()

val = [0.134202362865, 0.134202362865, 0.154144849556, 0.154144849556]
test.assertAlmostEqual(sfon[4, 1], val[rank])

test.printSummary()

FIN()
