# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
from code_aster.Utilities import PETSc

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI


CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()

MAIL = CA.ParallelMesh()
MAIL.readMedFile("mesh004b/%d.med" % rank, partitioned=True)
# MAIL.debugPrint()

MATER = DEFI_MATERIAU(THER=_F(LAMBDA=6.0e9, RHO_CP=1.0))

affectMat = CA.MaterialField(MAIL)
affectMat.addMaterialOnMesh(MATER)
affectMat.build()

MODT = AFFE_MODELE(
    MAILLAGE=MAIL,
    AFFE=_F(TOUT="OUI", PHENOMENE="THERMIQUE", MODELISATION="PLAN"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)


# charCine = CA.ThermalDirichletBC(MODT)
# charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Temp, 0., "EncastN")
# charCine.build()
charCine = AFFE_CHAR_CINE(MODELE=MODT, THER_IMPO=_F(TEMP=0, GROUP_NO="EncastN"))

CHT1 = AFFE_CHAR_THER(MODELE=MODT, FLUX_REP=_F(GROUP_MA="Press", FLUN=10.0), INFO=1)

# study = CA.PhysicalProblem(MODT, affectMat)
# study.addDirichletBC(charCine)
# study.addLoad(CHT1)
# dProblem = CA.DiscreteComputation(study)

vect_elem = CALC_VECT_ELEM(OPTION="CHAR_THER", CHARGE=CHT1)
matr_elem = CALC_MATR_ELEM(OPTION="RIGI_THER", MODELE=MODT, CHAM_MATER=affectMat)

numeDDL = CA.ParallelDOFNumbering()
numeDDL.computeNumbering([matr_elem])
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")
# numeDDL.debugPrint()

matrAsse = CA.AssemblyMatrixTemperatureReal()
matrAsse.setDOFNumbering(numeDDL)
matrAsse.assemble(matr_elem, charCine)
test.assertEqual(matrAsse.getType(), "MATR_ASSE_TEMP_R")
# matrAsse.debugPrint()

# retour = vect_elem.assembleWithLoadFunctions( numeDDL )
vecass = ASSE_VECTEUR(VECT_ELEM=vect_elem, NUME_DDL=numeDDL)
vcine = CALC_CHAR_CINE(NUME_DDL=numeDDL, CHAR_CINE=charCine)

matrAsse = FACTORISER(reuse=matrAsse, MATR_ASSE=matrAsse, METHODE="PETSC", PRE_COND="JACOBI")

retour = RESOUDRE(MATR=matrAsse, CHAM_NO=vecass, CHAM_CINE=vcine, ALGORITHME="CR", RESI_RELA=1e-9)

A = matrAsse.toPetsc()

v = PETSc.Viewer()
A.view(v)
v = PETSc.Viewer().createASCII("test.txt")
v.pushFormat(PETSc.Viewer.Format.ASCII_MATLAB)

rank = A.getComm().getRank()
print("rank=", rank)
rs, re = A.getOwnershipRange()
ce, _ = A.getSize()
rows = N.array(list(range(rs, re)), dtype=PETSc.IntType)
cols = N.array(list(range(0, ce)), dtype=PETSc.IntType)
rows = PETSc.IS().createGeneral(rows, comm=A.getComm())
cols = PETSc.IS().createGeneral(cols, comm=A.getComm())
(S,) = A.createSubMatrices(rows, cols)
v = PETSc.Viewer().createASCII("mesh004c_rank" + str(rank) + ".out", comm=S.getComm())
S.view(v)

# Use Dualized BC (with AFFE_CHAR_THER)
charTher = AFFE_CHAR_THER(MODELE=MODT, TEMP_IMPO=_F(TEMP=0, GROUP_NO="EncastN"))
vect_elem = CALC_VECT_ELEM(OPTION="CHAR_THER", CHARGE=CHT1)
matr_elem = CALC_MATR_ELEM(OPTION="RIGI_THER", CHARGE=charTher, MODELE=MODT, CHAM_MATER=affectMat)

numeDDL = CA.ParallelDOFNumbering()
numeDDL.computeNumbering([matr_elem])
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")

matrAsse = CA.AssemblyMatrixTemperatureReal()
matrAsse.setDOFNumbering(numeDDL)
matrAsse.assemble(matr_elem)
test.assertEqual(matrAsse.getType(), "MATR_ASSE_TEMP_R")

vecass = ASSE_VECTEUR(VECT_ELEM=vect_elem, NUME_DDL=numeDDL)

matrAsse = FACTORISER(reuse=matrAsse, MATR_ASSE=matrAsse, METHODE="PETSC", PRE_COND="JACOBI")

retour = RESOUDRE(MATR=matrAsse, CHAM_NO=vecass, ALGORITHME="GCR", RESI_RELA=1e-9)

# Export / import to PETSc
U = retour

pU = U.toPetsc()
V = U.copy()

V.setValues(0.0)
test.assertEqual(V.norm("NORM_2"), 0)

test.assertEqual((U - V.fromPetsc(pU)).norm("NORM_2"), 0)

scaling = 1000.0
V.fromPetsc(pU, scaling)
for lag in numeDDL.getLagrangeDOFs(local=True):
    test.assertEqual((V[lag] - U[lag] * 1000.0), 0)

U.applyLagrangeScaling(scaling)
test.assertEqual((U - V).norm("NORM_2"), 0)

test.printSummary()

FIN()
