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

import numpy as N
import petsc4py
import code_aster
from code_aster.Commands import *
from code_aster import MPI


code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

MAIL = code_aster.ParallelMesh()
MAIL.readMedFile("mesh004b/%d.med" % rank, partitioned=True)
# MAIL.debugPrint()

MATER = DEFI_MATERIAU(THER=_F(LAMBDA=6.0e9, RHO_CP=1.0))

affectMat = code_aster.MaterialField(MAIL)
affectMat.addMaterialOnMesh(MATER)
affectMat.build()

MODT = AFFE_MODELE(
    MAILLAGE=MAIL,
    AFFE=_F(TOUT="OUI", PHENOMENE="THERMIQUE", MODELISATION="PLAN"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)


# charCine = code_aster.ThermalDirichletBC(MODT)
# charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Temp, 0., "EncastN")
# charCine.build()
charCine = AFFE_CHAR_CINE(MODELE=MODT, THER_IMPO=_F(TEMP=0, GROUP_NO="EncastN"))

CHT1 = AFFE_CHAR_THER(MODELE=MODT, FLUX_REP=_F(GROUP_MA="Press", FLUN=10.0), INFO=1)

# study = code_aster.PhysicalProblem(MODT, affectMat)
# study.addDirichletBC(charCine)
# study.addLoad(CHT1)
# dProblem = code_aster.DiscreteComputation(study)

vect_elem = CALC_VECT_ELEM(OPTION="CHAR_THER", CHARGE=CHT1)
matr_elem = CALC_MATR_ELEM(OPTION="RIGI_THER", MODELE=MODT, CHAM_MATER=affectMat)

numeDDL = code_aster.ParallelDOFNumbering()
numeDDL.setElementaryMatrix(matr_elem)
numeDDL.computeNumbering()
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")
# numeDDL.debugPrint()

matrAsse = code_aster.AssemblyMatrixTemperatureReal()
matrAsse.addElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.addDirichletBC(charCine)
matrAsse.assemble()
test.assertEqual(matrAsse.getType(), "MATR_ASSE_TEMP_R")
# matrAsse.debugPrint()

# retour = vect_elem.assembleWithLoadFunctions( numeDDL )
vecass = ASSE_VECTEUR(VECT_ELEM=vect_elem, NUME_DDL=numeDDL)
vcine = CALC_CHAR_CINE(NUME_DDL=numeDDL, CHAR_CINE=charCine)

matrAsse = FACTORISER(reuse=matrAsse, MATR_ASSE=matrAsse, METHODE="PETSC", PRE_COND="JACOBI")

retour = RESOUDRE(MATR=matrAsse, CHAM_NO=vecass, CHAM_CINE=vcine, ALGORITHME="CR", RESI_RELA=1e-9)

A = matrAsse.toPetsc()

v = petsc4py.PETSc.Viewer()
A.view(v)
v = petsc4py.PETSc.Viewer().createASCII("test.txt")
v.pushFormat(petsc4py.PETSc.Viewer.Format.ASCII_MATLAB)

rank = A.getComm().getRank()
print("rank=", rank)
rs, re = A.getOwnershipRange()
ce, _ = A.getSize()
rows = N.array(list(range(rs, re)), dtype=petsc4py.PETSc.IntType)
cols = N.array(list(range(0, ce)), dtype=petsc4py.PETSc.IntType)
rows = petsc4py.PETSc.IS().createGeneral(rows, comm=A.getComm())
cols = petsc4py.PETSc.IS().createGeneral(cols, comm=A.getComm())
(S,) = A.createSubMatrices(rows, cols)
v = petsc4py.PETSc.Viewer().createASCII("mesh004c_rank" + str(rank) + ".out", comm=S.getComm())
S.view(v)

# Use Dualized BC (with AFFE_CHAR_THER)
charTher = AFFE_CHAR_THER(MODELE=MODT, TEMP_IMPO=_F(TEMP=0, GROUP_NO="EncastN"))
vect_elem = CALC_VECT_ELEM(OPTION="CHAR_THER", CHARGE=CHT1)
matr_elem = CALC_MATR_ELEM(OPTION="RIGI_THER", CHARGE=charTher, MODELE=MODT, CHAM_MATER=affectMat)

numeDDL = code_aster.ParallelDOFNumbering()
numeDDL.setElementaryMatrix(matr_elem)
numeDDL.computeNumbering()
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")

matrAsse = code_aster.AssemblyMatrixTemperatureReal()
matrAsse.addElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.assemble()
test.assertEqual(matrAsse.getType(), "MATR_ASSE_TEMP_R")

vecass = ASSE_VECTEUR(VECT_ELEM=vect_elem, NUME_DDL=numeDDL)

matrAsse = FACTORISER(reuse=matrAsse, MATR_ASSE=matrAsse, METHODE="PETSC", PRE_COND="JACOBI")

retour = RESOUDRE(MATR=matrAsse, CHAM_NO=vecass, ALGORITHME="GCR", RESI_RELA=1e-9)

# Export / import to PETSc
test = code_aster.TestCase()
U = retour

pU = U.toPetsc(numeDDL)
V = U.duplicate()

V.setValues(0.0)
test.assertEqual(V.norm("NORM_2"), 0)

V.fromPetsc(numeDDL, pU)
test.assertEqual((U - V).norm("NORM_2"), 0)

scaling = 1000.0
V.fromPetsc(numeDDL, pU, scaling)
for lag in numeDDL.getRowsAssociatedToLagrangeMultipliers(local=True):
    test.assertEqual((V[lag] - U[lag] * 1000.0), 0)

U.applyLagrangeScaling(scaling)
test.assertEqual((U - V).norm("NORM_2"), 0)

test.printSummary()

FIN()
