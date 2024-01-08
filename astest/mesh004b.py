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

import numpy as N
from code_aster import CA
from code_aster.Commands import *
from code_aster.CA import MPI
from code_aster.Utilities import PETSc

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

pMesh = CA.ParallelMesh()
pMesh.readMedFile("mesh004b/%d.med" % rank, partitioned=True)

MATER = DEFI_MATERIAU(ELAS=_F(E=10000.0, NU=0.0, RHO=1.0))

affectMat = CA.MaterialField(pMesh)
affectMat.addMaterialOnMesh(MATER)
affectMat.build()

MODT = AFFE_MODELE(
    MAILLAGE=pMesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN")
)

# MODT = CA.Model(MAIL))

charCine = CA.MechanicalDirichletBC(MODT)
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dx, 0.0, "EncastN")
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dy, 0.0, "EncastN")
charCine.build()

# CHT1 = AFFE_CHAR_MECA(MODELE=MODT,
#                      PESANTEUR=_F(GRAVITE=1.0,
#                                   DIRECTION=(0.0, -1.0, 0.0),),
#                      INFO=1,
#                      VERI_NORM='NON',)


CHT1 = AFFE_CHAR_MECA(MODELE=MODT, PRES_REP=_F(GROUP_MA="Press", PRES=-10), INFO=1, VERI_NORM="NON")
vect_elem = CALC_VECT_ELEM(OPTION="CHAR_MECA", CHARGE=CHT1)

study = CA.PhysicalProblem(MODT, affectMat)
study.addDirichletBC(charCine)
study.addLoad(CHT1)
study.computeDOFNumbering()
dComputation = CA.DiscreteComputation(study)
matr_elem = dComputation.getElasticStiffnessMatrix()

monSolver = CA.PetscSolver(RENUM="SANS", PRE_COND="SANS")

numeDDL = CA.ParallelDOFNumbering()
numeDDL.computeNumbering(MODT, study.getListOfLoads())
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")

# compute Neumman
vecass = ASSE_VECTEUR(VECT_ELEM=vect_elem, NUME_DDL=numeDDL)
print("vecass=", vecass.getValues())
retour = dComputation.getNeumannForces(0, 0, 0)


matrAsse = CA.AssemblyMatrixDisplacementReal()
matrAsse.addElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.addDirichletBC(charCine)
matrAsse.assemble()
test.assertEqual(matrAsse.getType(), "MATR_ASSE_DEPL_R")

print("retour=", retour.getValues())

monSolver.factorize(matrAsse)
resu = monSolver.solve(retour)


# ------------------------------------
# tests in local numbering
physicalRows = numeDDL.getPhysicalDOFs(local=True)
test.assertListEqual(
    physicalRows, [i * 2 + j for i in range(pMesh.getNumberOfNodes()) for j in range(2)]
)
multipliersRows = numeDDL.getLagrangeDOFs(local=True)
test.assertListEqual(multipliersRows, [])
test.assertFalse(numeDDL.useLagrangeDOF())
test.assertFalse(numeDDL.useSingleLagrangeDOF())
test.assertEqual(numeDDL.getComponents(), ["DX", "DY"])
test.assertEqual(numeDDL.getComponentFromNode(0, local=True), ["DX", "DY"])
test.assertEqual(numeDDL.getNodeFromDOF(0, local=True), 0)
test.assertTrue(numeDDL.isPhysicalDOF(0, local=True))
test.assertEqual(numeDDL.getNumberOfDOFs(local=True), 2 * pMesh.getNumberOfNodes())
test.assertEqual(numeDDL.getNumberOfDOFs(local=False), 16)
test.assertEqual(numeDDL.getPhysicalQuantity(), "DEPL_R")
ghostRows = numeDDL.getGhostDOFs(local=True)
test.assertListEqual(ghostRows, [[6, 7, 10, 11], [4, 5, 8, 9]][rank])


# ------------------------------------
# tests in global numbering
physicalRows = numeDDL.getPhysicalDOFs(local=False)
test.assertListEqual(
    physicalRows,
    [
        numeDDL.localToGlobalDOF(d)
        for d in [i * 2 + j for i in range(pMesh.getNumberOfNodes()) for j in range(2)]
    ],
)
test.assertListEqual(
    physicalRows,
    [[0, 1, 2, 3, 4, 5, 12, 13, 6, 7, 14, 15], [8, 9, 10, 11, 4, 5, 12, 13, 6, 7, 14, 15]][rank],
)

ghostRows = numeDDL.getGhostDOFs(local=False)
test.assertListEqual(
    ghostRows, [numeDDL.localToGlobalDOF(i) for i in [[6, 7, 10, 11], [4, 5, 8, 9]][rank]]
)
test.assertListEqual(ghostRows, [[12, 13, 14, 15], [4, 5, 6, 7]][rank])


# ------------------------------------
# Some petsc4py manipulations
pA = matrAsse.toPetsc()
v = PETSc.Viewer().createASCII("mesh004b.out")
v.pushFormat(PETSc.Viewer.Format.ASCII_DENSE)
pA.view(v)

rank = pA.getComm().getRank()
print("rank=", rank)
rs, re = pA.getOwnershipRange()
ce, _ = pA.getSize()
rows = N.array(list(range(rs, re)), dtype=PETSc.IntType)
cols = N.array(list(range(0, ce)), dtype=PETSc.IntType)
rows = PETSc.IS().createGeneral(rows, comm=pA.getComm())
cols = PETSc.IS().createGeneral(cols, comm=pA.getComm())
(S,) = pA.createSubMatrices(rows, cols)
v = PETSc.Viewer().createASCII("mesh004b_rank" + str(rank) + ".out", comm=S.getComm())
S.view(v)


# ------------------------------------
# Scaling validation
from code_aster.LinearAlgebra import MatrixScaler
from code_aster.Utilities import logger

pA_unscaled = matrAsse.toPetsc()
pA_unscaled.view()
# The Scling object
S = MatrixScaler.MatrixScaler()
logger.setLevel(2)
newMat = matrAsse.copy()
# Compute scaling with DX and DY gathered (default behavior)
S.computeScaling(matrAsse)
S.scaleMatrix(newMat)

pA_scaled = newMat.toPetsc()
pA_scaled.view()
nt = PETSc.NormType.NORM_INFINITY
test.assertAlmostEqual(pA_unscaled.norm(nt), 43055.55555560758)
test.assertAlmostEqual(pA_scaled.norm(nt), 1.0)


newMat = matrAsse.copy()
# Compute scaling with DX and DY separately
S.computeScaling(matrAsse, merge_dof=[])
S.scaleMatrix(newMat)

pA_scaled = newMat.toPetsc()
pA_scaled.view()
nt = PETSc.NormType.NORM_INFINITY
test.assertAlmostEqual(pA_unscaled.norm(nt), 43055.55555560758)
test.assertAlmostEqual(pA_scaled.norm(nt), 1.0)

rhs = vecass.copy()
init_norm = rhs.norm("NORM_INFINITY")
S.scaleRHS(rhs)
test.assertAlmostEqual(rhs.norm("NORM_INFINITY"), 286.21852876537895)


sol = vecass.copy()
S.unscaleSolution(rhs)
test.assertAlmostEqual(rhs.norm("NORM_INFINITY"), init_norm)


logger.setLevel(0)

test.printSummary()

FIN()
