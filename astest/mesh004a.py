# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
from code_aster.Utilities import petscInitialize
from code_aster import MPI


code_aster.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

# force PETSc to start before solves for testing purpose only - no need in regular study
petscInitialize("-ksp_view -log_view -ksp_monitor")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
print("Nb procs", MPI.ASTER_COMM_WORLD.Get_size())
print("Rank", MPI.ASTER_COMM_WORLD.Get_rank())

pMesh2 = code_aster.ParallelMesh()
pMesh2.readMedFile("mesh004a/%d.med" % rank, partitioned=True)
pMesh2 = DEFI_GROUP(reuse=pMesh2, MAILLAGE=pMesh2, CREA_GROUP_NO=_F(TOUT_GROUP_MA="OUI"))
del pMesh2

pMesh = code_aster.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med" % rank, partitioned=True)
pMesh.debugPrint(rank + 30)

model = code_aster.Model(pMesh)
test.assertEqual(model.getType(), "MODELE_SDASTER")
model.addModelingOnMesh(code_aster.Physics.Mechanics, code_aster.Modelings.Tridimensional)
model.build()

testMesh = model.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

model.debugPrint(rank + 30)

acier = DEFI_MATERIAU(ELAS=_F(E=2.0e11, NU=0.3))
acier.debugPrint(8)

affectMat = code_aster.MaterialField(pMesh)

testMesh2 = affectMat.getMesh()
test.assertEqual(testMesh2.getType(), "MAILLAGE_P")

affectMat.addMaterialOnMesh(acier)
affectMat.build()

charCine = code_aster.MechanicalDirichletBC(model)
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dx, 0.0, "COTE_B")
charCine.build()

study = code_aster.PhysicalProblem(model, affectMat)
study.addDirichletBC(charCine)
study.computeDOFNumbering()
dComputation = code_aster.DiscreteComputation(study)
matr_elem = dComputation.getElasticStiffnessMatrix()

listLoads = study.getListOfLoads()

numeDDL = code_aster.ParallelDOFNumbering()
numeDDL.computeNumbering([matr_elem])
numeDDL.debugPrint(rank + 30)

matrAsse = code_aster.AssemblyMatrixDisplacementReal()
matrAsse.addElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.addDirichletBC(charCine)
matrAsse.assemble()
matrAsse.debugPrint(rank + 30)

ccid = matrAsse.getDirichletBCDOFs()
bcNb = {0: 142, 1: 148, 2: 0, 3: 0}
test.assertEqual(sum(ccid), bcNb[rank])
test.assertEqual(len(ccid), numeDDL.getNumberOfDOFs(local=True) + 1)

vec = code_aster.FieldOnNodesReal(model)
vec.setValues(1.0)
study.zeroDirichletBCDOFs(vec)
test.assertEqual(sum(vec.getValues()), numeDDL.getNumberOfDOFs(local=True) - sum(ccid[:-1]))

retour = matrAsse.getDOFNumbering()
test.assertEqual(retour.isParallel(), True)

# tests on numbering of DOFs
physicalRows = numeDDL.getPhysicalDOFs(local=True)
test.assertListEqual(physicalRows, list(range(3 * len(pMesh.getNodes(localNumbering=True)))))
multipliersRows = numeDDL.getLagrangeDOFs(local=True)
test.assertListEqual(multipliersRows, [])
test.assertFalse(numeDDL.useLagrangeDOF())
test.assertFalse(numeDDL.useSingleLagrangeDOF())
test.assertEqual(numeDDL.getComponents(), ["DX", "DY", "DZ"])
test.assertEqual(numeDDL.getComponentFromNode(0, local=True), ["DX", "DY", "DZ"])
test.assertEqual(numeDDL.getNodeFromDOF(0, local=True), 0)
test.assertTrue(numeDDL.isPhysicalDOF(0, local=True))
test.assertEqual(numeDDL.getNumberOfDOFs(local=True), 3 * len(pMesh.getNodes(localNumbering=True)))
test.assertEqual(numeDDL.getNumberOfDOFs(local=False), 3993)
test.assertEqual(numeDDL.getPhysicalQuantity(), "DEPL_R")
nnodes = pMesh.getNumberOfNodes()
test.assertEqual(
    numeDDL.getNodeAndComponentFromDOF(local=True)[::3], [(i, "DX") for i in range(nnodes)]
)
l2G = pMesh.getLocalToGlobalNodeIds()
test.assertEqual(
    numeDDL.getNodeAndComponentFromDOF(local=False)[::3], [(l2G[i], "DX") for i in range(nnodes)]
)
test.printSummary()

code_aster.close()
