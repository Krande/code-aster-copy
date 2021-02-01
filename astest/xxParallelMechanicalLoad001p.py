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

# force PETSc to start before solves for testing purpose only - no need in regular study
LinearAlgebra.petscInitialize("-ksp_view -log_view -ksp_monitor")

test = code_aster.TestCase()

code_aster.init("--test")

rank = code_aster.getMPIRank()

pMesh2 = code_aster.ParallelMesh()
pMesh2.readMedFile("xxParallelMechanicalLoad001b/%d.med"%rank, True)

model = AFFE_MODELE(MAILLAGE = pMesh2,
                    AFFE = _F(MODELISATION = "D_PLAN_INCO_UPG",
                              PHENOMENE = "MECANIQUE",
                              TOUT = "OUI",),)

char_cin = AFFE_CHAR_CINE(MODELE=model,
                          MECA_IMPO=(_F(GROUP_NO="N2",
                                        DX=0.,DY=0.,DZ=0.,),
                                     _F(GROUP_NO="N4",
                                        DX=0.,DY=0.,DZ=0.,),),)

char_meca = AFFE_CHAR_MECA(MODELE=model,
                           LIAISON_DDL=_F(GROUP_NO=("N1", "N3"),
                                          DDL=('PRES','PRES'),
                                          COEF_MULT=(1.0,1.0),
                                          COEF_IMPO=0,),
                           DDL_IMPO=(_F(GROUP_NO="N1",PRES=200000.0),
                                          ))

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0,
                               NU=0.499999,),)

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2,
                       AFFE=_F(TOUT='OUI',
                               MATER=MATER1,),)

# =======================================================================


study = code_aster.StudyDescription(model, AFFMAT)
study.addDirichletBC(char_cin)
study.addParallelMechanicalLoad(char_meca)
dProblem = code_aster.DiscreteProblem(study)
vect_elem = dProblem.buildElementaryMechanicalLoadsVector()
matr_elem = dProblem.computeMechanicalStiffnessMatrix()

numeDDL = code_aster.ParallelDOFNumbering()
numeDDL.setElementaryMatrix(matr_elem)
numeDDL.computeNumbering()
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")

matrAsse = code_aster.AssemblyMatrixDisplacementReal()
matrAsse.appendElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.addDirichletBC(char_cin)
matrAsse.build()

matrAsse.print()

test.assertRaises(RuntimeError, lambda: list(numeDDL.getRowsAssociatedToPhysicalDofs()))
test.assertRaises(RuntimeError, lambda: list(numeDDL.getRowsAssociatedToLagrangeMultipliers()))
test.assertRaises(RuntimeError, lambda: list(numeDDL.getComponentsAssociatedToNode(1)))
test.assertRaises(RuntimeError, lambda: list(numeDDL.getNodeAssociatedToRow(1)))

physicalRows = numeDDL.getRowsAssociatedToPhysicalDofs(local=True)
if rank == 0: test.assertListEqual(physicalRows, [1,2,3,4,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40])
if rank == 1: test.assertListEqual(physicalRows, [1,2,3,4,7,8,9,10,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42])

multipliersRows = numeDDL.getRowsAssociatedToLagrangeMultipliers(local=True)
if rank == 0: test.assertListEqual(multipliersRows, [5, 10])
if rank == 1: test.assertListEqual(multipliersRows, [5, 6, 11, 12])

test.assertTrue(numeDDL.useLagrangeMultipliers())
test.assertFalse(numeDDL.useSingleLagrangeMultipliers())
test.assertEqual(numeDDL.getComponents(), ['DX', 'DY','GONF','LAGR', 'PRES'])
test.assertEqual(numeDDL.getComponentsAssociatedToNode(1, local=True), ['DX', 'DY', 'PRES', 'GONF'])
test.assertEqual(numeDDL.getNodeAssociatedToRow(1, local=True), 1)
if rank == 0: test.assertEqual(numeDDL.getNodeAssociatedToRow(5, local=True), 0)
if rank == 1: test.assertEqual(numeDDL.getNodeAssociatedToRow(5, local=True), -2)
if rank == 0: test.assertEqual(numeDDL.getComponentAssociatedToRow(5, local=True), ' ')
if rank == 1: test.assertEqual(numeDDL.getComponentAssociatedToRow(5, local=True), 'PRES')
if rank == 0: test.assertEqual(numeDDL.getNumberOfDofs(local=True), 40)
if rank == 1: test.assertEqual(numeDDL.getNumberOfDofs(local=True), 42)
test.assertEqual(numeDDL.getNumberOfDofs(local=False), 56)
test.assertEqual(numeDDL.getPhysicalQuantity(), 'DEPL_R')


# Solve

monSolver = code_aster.PetscSolver(code_aster.Renumbering.Sans)
monSolver.setPreconditioning(code_aster.Preconditioning.SimplePrecisionLdlt)
monSolver.setLowRankThreshold(1.e-6)
monSolver.setSolverResidual(1.e-10)

monSolver.matrixFactorization(matrAsse)
test.assertEqual(matrAsse.getType(), "MATR_ASSE_DEPL_R")

vcine = dProblem.buildKinematicsLoad(numeDDL, 0.)
test.assertEqual(vect_elem.getType(), "VECT_ELEM_DEPL_R")
vasse = vect_elem.assembleVector(numeDDL)
resu = monSolver.solveRealLinearSystemWithKinematicsLoad(matrAsse, vcine, vasse)

sfon = resu.exportToSimpleFieldOnNodes()
sfon.updateValuePointers()


# DX displacement on nodes "N1" and "N3", comparison with sequential results
if rank == 0:
    test.assertAlmostEqual(sfon.getValue(1, 0), 1.14977255749554, 6)
elif rank == 1:
    test.assertAlmostEqual(sfon.getValue(1, 0), 1.14977255749554, 6)


test.assertRaises(RuntimeError, lambda: list(vasse.setDirichletBC(GROUP_NO=('N1'),PRES=1000,DX=10.)))


test.printSummary()


code_aster.close()
