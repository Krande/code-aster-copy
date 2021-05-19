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
from code_aster import LinearAlgebra
from code_aster import MPI


# force PETSc to start before solves for testing purpose only - no need in regular study
LinearAlgebra.petscInitialize("-ksp_view -log_view -ksp_monitor")

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.COMM_WORLD.Get_rank()
print("Nb procs", MPI.COMM_WORLD.Get_size())
print("Rank", MPI.COMM_WORLD.Get_rank())

pMesh2 = code_aster.ParallelMesh()
pMesh2.readMedFile("mesh004a/%d.med"%rank, True)
pMesh2=DEFI_GROUP(reuse =pMesh2,MAILLAGE=pMesh2,CREA_GROUP_NO=_F(TOUT_GROUP_MA='OUI'))
del pMesh2

pMesh = code_aster.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med"%rank, True)
pMesh.debugPrint(rank+30)

model = code_aster.Model(pMesh)
test.assertEqual(model.getType(), "MODELE_SDASTER")
model.addModelingOnMesh(code_aster.Physics.Mechanics,
                           code_aster.Modelings.Tridimensional)
model.build()

testMesh = model.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

model.debugPrint(rank+30)

acier = DEFI_MATERIAU(ELAS = _F(E = 2.e11,
                                NU = 0.3,),)
acier.debugPrint(8)

affectMat = code_aster.MaterialField(pMesh)

testMesh2 = affectMat.getMesh()
test.assertEqual(testMesh2.getType(), "MAILLAGE_P")

affectMat.addMaterialsOnMesh(acier)
affectMat.buildWithoutExternalVariable()

charCine = code_aster.MechanicalDirichletBC(model)
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dx, 0., "COTE_B")
charCine.build()

study = code_aster.StudyDescription(model, affectMat)
dProblem = code_aster.DiscreteProblem(study)
matr_elem = dProblem.computeMechanicalStiffnessMatrix()

monSolver = code_aster.MumpsSolver(code_aster.Renumbering.Metis)

numeDDL = code_aster.ParallelDOFNumbering()
numeDDL.setElementaryMatrix(matr_elem)
numeDDL.computeNumbering()
numeDDL.debugPrint(rank+30)

matrAsse = code_aster.AssemblyMatrixDisplacementReal()
matrAsse.appendElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.addDirichletBC(charCine)
matrAsse.build()
matrAsse.debugPrint(rank+30)

retour = matrAsse.getDOFNumbering()
test.assertEqual(retour.isParallel(), True)

test.assertRaises(RuntimeError, lambda: list(numeDDL.getRowsAssociatedToPhysicalDofs()))
test.assertRaises(RuntimeError, lambda: list(numeDDL.getRowsAssociatedToLagrangeMultipliers()))
test.assertRaises(RuntimeError, lambda: list(numeDDL.getComponentsAssociatedToNode(1)))
test.assertRaises(RuntimeError, lambda: list(numeDDL.getNodeAssociatedToRow(1)))

physicalRows = numeDDL.getRowsAssociatedToPhysicalDofs(local=True)
test.assertListEqual(physicalRows, list(range(1,3*len(pMesh.getNodes(localNumbering=True))+1)))
multipliersRows = numeDDL.getRowsAssociatedToLagrangeMultipliers(local=True)
test.assertListEqual(multipliersRows, [])
test.assertFalse(numeDDL.useLagrangeMultipliers())
test.assertFalse(numeDDL.useSingleLagrangeMultipliers())
test.assertEqual(numeDDL.getComponents(), ['DX', 'DY', 'DZ'])
test.assertEqual(numeDDL.getComponentsAssociatedToNode(1, local=True), ['DX', 'DY', 'DZ'])
test.assertEqual(numeDDL.getNodeAssociatedToRow(1, local=True), 1)
test.assertEqual(numeDDL.getNumberOfDofs(local=True), 3*len(pMesh.getNodes(localNumbering=True)))
test.assertEqual(numeDDL.getNumberOfDofs(local=False), 3993)
test.assertEqual(numeDDL.getPhysicalQuantity(), 'DEPL_R')

test.printSummary()

code_aster.close()
