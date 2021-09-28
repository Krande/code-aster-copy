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
from code_aster import MPI


code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.COMM_WORLD.Get_rank()
pMesh = code_aster.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med"%rank, True)
DEFI_GROUP(reuse=pMesh,
           MAILLAGE=pMesh,
           CREA_GROUP_MA=(_F(NOM='BLABLA',
                             OPTION='SPHERE',
                             POINT=(0.2, 0.2, 0.2),
                             RAYON = 0.2),),)
pMesh.debugPrint(30+rank)


DEFI_GROUP(reuse=pMesh,
           MAILLAGE=pMesh,
           CREA_GROUP_MA=(_F(NOM='BLABLA',
                             OPTION='SPHERE',
                             POINT=(0.2, 0.2, 0.2),
                             RAYON = 0.2),),)

list_cells = pMesh.getCells( 'BLABLA' )
nb_cells = [ 72 , 4 , 0, 4]
test.assertEqual(len(list_cells), nb_cells[rank])

monModel = code_aster.Model(pMesh)
monModel.addModelingOnMesh(code_aster.Physics.Mechanics,
                              code_aster.Modelings.Tridimensional)
monModel.build()

testMesh = monModel.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

acier = DEFI_MATERIAU(ELAS = _F(E = 2.e11,
                                NU = 0.3,),)

affectMat = code_aster.MaterialField(pMesh)
affectMat.addMaterialsOnMesh( acier )
affectMat.buildWithoutExternalStateVariables()

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

monSolver = code_aster.PetscSolver(code_aster.Renumbering.Sans)
monSolver.setPreconditioning(code_aster.Preconditioning.Sor)

mecaStatique = code_aster.LinearStaticAnalysis(monModel, affectMat)
mecaStatique.addDirichletBC(charCine)
mecaStatique.addDirichletBC(charCine2)
mecaStatique.setLinearSolver(monSolver)

resu = mecaStatique.execute()

resu=CALC_CHAMP(RESULTAT=resu, reuse=resu, CONTRAINTE=('SIEF_ELGA'))

IMPR_RESU(FICHIER_UNIQUE='OUI',
          FORMAT='MED',
          RESU=_F(RESULTAT=resu,),
          VERSION_MED='4.0.0')

IMPR_RESU(FICHIER_UNIQUE='OUI',
          FORMAT='MED',
          RESU=_F(RESULTAT=resu, GROUP_NO="COTE_H"),
          VERSION_MED='4.0.0',
          UNITE=81)

IMPR_RESU(FICHIER_UNIQUE='OUI',
          FORMAT='MED',
          RESU=_F(RESULTAT=resu, GROUP_MA="BLABLA"),
          VERSION_MED='4.0.0',
          UNITE=82)

#resu.printMedFile("fort."+str(rank+40)+".med")

MyFieldOnNodes = resu.getFieldOnNodesReal("DEPL", 0)
sfon = MyFieldOnNodes.exportToSimpleFieldOnNodes()
sfon.updateValuePointers()

val = [0.134228076192 , 0.134176297047, 0.154099687654, 0.154189676715]
test.assertAlmostEqual(sfon.getValue(4, 1), val[rank])

test.printSummary()

FIN()
