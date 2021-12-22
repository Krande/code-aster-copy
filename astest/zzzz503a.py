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
import libaster
import numpy as np

code_aster.init("--test")

test = code_aster.TestCase()

monMaillage = code_aster.Mesh()
monMaillage.readMedFile("zzzz503a.mmed")

monModel = code_aster.Model(monMaillage)
monModel.addModelingOnMesh(code_aster.Physics.Mechanics, code_aster.Modelings.Tridimensional)
monModel.build()
test.assertEqual(monModel.getType(), "MODELE_SDASTER")

YOUNG = 200000.0
POISSON = 0.3

Kinv = 3.2841e-4
Kv = 1. / Kinv
SY = 437.0
Rinf = 758.0
Qzer = 758.0 - 437.
Qinf = Qzer + 100.
b = 2.3
C1inf = 63767.0 / 2.0
C2inf = 63767.0 / 2.0
Gam1 = 341.0
Gam2 = 341.0
C_Pa = 1.e+6

acier = DEFI_MATERIAU(ELAS=_F(E=YOUNG,
                              NU=POISSON, ),
                      VISCOCHAB=_F(K=SY * C_Pa,
                                   B=b,
                                   MU=10,
                                   Q_M=Qinf * C_Pa,
                                   Q_0=Qzer * C_Pa,
                                   C1=C1inf * C_Pa,
                                   C2=C2inf * C_Pa,
                                   G1_0=Gam1,
                                   G2_0=Gam2,
                                   K_0=Kv * C_Pa,
                                   N=11,
                                   A_K=1., ), )
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


study = code_aster.PhysicalProblem(monModel, affectMat)
study.addLoad(CharMeca1)
study.addLoad(CharMeca2)
listLoads = study.getListOfLoads()
study.computeDOFNumbering()
dComputation = code_aster.DiscreteComputation(study)
# compute Neumann
retour = dComputation.neumann([1, 0, 0])
matr_elem = dComputation.computeMechanicalStiffnessMatrix()

test.assertEqual(matr_elem.getType(), "MATR_ELEM_DEPL_R")

monSolver = code_aster.MumpsSolver()

numeDDL = code_aster.DOFNumbering()
numeDDL.setElementaryMatrix(matr_elem)
numeDDL.computeNumbering()
print(numeDDL.getListOfLoads(), flush=True)
study.setDOFNumbering(numeDDL)
# numeDDL.debugPrint(6)
test.assertEqual(numeDDL.getType(), "NUME_DDL_SDASTER")
test.assertFalse(numeDDL.hasDirichletBC())

ccid = numeDDL.getDirichletBCDOFs()
test.assertEqual(sum(ccid), 0)
test.assertEqual(len(ccid), numeDDL.getNumberOfDofs())

matrAsse = code_aster.AssemblyMatrixDisplacementReal()
matrAsse.appendElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.build()

x = matrAsse.EXTR_MATR(sparse=True)
test.assertTrue('numpy' in str(type(x[0])))

# test setValues
# -----------------
values, idx, jdx, neq = matrAsse.EXTR_MATR(sparse=True)
K1 = matrAsse.EXTR_MATR()

neq = K1.shape[0]
matrAsse.setValues([0, 1], [0, 1], [1., 1.])
K2 = matrAsse.EXTR_MATR()
test.assertAlmostEqual(np.linalg.norm(K2), np.sqrt(2))

matrAsse.setValues(idx.tolist(), jdx.tolist(), values.tolist())
K3 = matrAsse.EXTR_MATR()
test.assertEqual(np.linalg.norm(K1 - K3), 0)

try:
    import petsc4py

    A = matrAsse.toPetsc()
except (ImportError, NotImplementedError):
    pass
else:
    v = petsc4py.PETSc.Viewer()
    A.view(v)
    v = petsc4py.PETSc.Viewer().createASCII("test.txt")
    v.pushFormat(petsc4py.PETSc.Viewer.Format.ASCII_MATLAB)
    A.view(v)

monSolver.factorize(matrAsse)
test.assertEqual(matrAsse.getType(), "MATR_ASSE_DEPL_R")

vcine = dComputation.dirichletBC(0.)
resu = monSolver.solve(retour, vcine)

y = resu.EXTR_COMP()
test.assertEqual(len(y.valeurs), 81)

resu2 = resu.exportToSimpleFieldOnNodes()
resu2.updateValuePointers()
test.assertAlmostEqual(resu2.getValue(6, 0), 0.000757555469653289)

resu.printMedFile("fort.med")

# test setValues + solve
# ----------------------
matrAsse.setValues(idx.tolist(), jdx.tolist(), [10 * v for v in values])
monSolver.factorize(matrAsse)
resu = monSolver.solve(retour)
resu2 = resu.exportToSimpleFieldOnNodes()
resu2.updateValuePointers()
test.assertAlmostEqual(resu2.getValue(6, 0), 0.000757555469653289 / 10.)

# To be sure that vcine is Permanent #30689
libaster.deleteTemporaryObjects()
test.assertTrue(vcine.updateValuePointers())

# check other solvers attributes
test.assertEqual(monSolver.getSolverName(), "MUMPS")
test.assertTrue(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = code_aster.LdltSolver()
test.assertEqual(monSolver.getSolverName(), "LDLT")
test.assertFalse(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = code_aster.MultFrontSolver()
test.assertEqual(monSolver.getSolverName(), "MULT_FRONT")
test.assertFalse(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = code_aster.PetscSolver()
test.assertEqual(monSolver.getSolverName(), "PETSC")
test.assertTrue(monSolver.supportParallelMesh(), "support of ParallelMesh")

monSolver = code_aster.GcpcSolver()
test.assertEqual(monSolver.getSolverName(), "GCPC")
test.assertFalse(monSolver.supportParallelMesh(), "support of ParallelMesh")


test.printSummary()

FIN()
