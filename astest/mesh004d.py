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

import numpy as N
import code_aster
from code_aster.Commands import *

code_aster.init("--test")

test = code_aster.TestCase()

rank = code_aster.getMPIRank()

MAIL = code_aster.ParallelMesh()
MAIL.readMedFile("mesh004b/%d.med" % rank, True)
# MAIL.debugPrint()

MATER = DEFI_MATERIAU(THER=_F(LAMBDA=6.E9, RHO_CP=1.))

affectMat = code_aster.MaterialField(MAIL)
affectMat.addMaterialsOnMesh(MATER)
affectMat.buildWithoutExternalVariable()

MODT = AFFE_MODELE(MAILLAGE=MAIL,
                   AFFE=_F(TOUT='OUI',
                           PHENOMENE='THERMIQUE',
                           MODELISATION='PLAN',),
                   DISTRIBUTION=_F(METHODE='CENTRALISE',),)


charCine = code_aster.ThermalDirichletBC(MODT)
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Temp, 0., "EncastN")
charCine.build()
charCine.debugPrint()

CHT1 = AFFE_CHAR_THER(MODELE=MODT,
                      FLUX_REP=_F(GROUP_MA='Press', FLUN=10.,),
                      INFO=1,)

# study = code_aster.StudyDescription(MODT, affectMat)
# study.addDirichletBC(charCine)
# study.addLoad(CHT1)
# dProblem = code_aster.DiscreteProblem(study)

vect_elem = CALC_VECT_ELEM(OPTION='CHAR_THER', CHARGE=CHT1)
matr_elem = CALC_MATR_ELEM(OPTION='RIGI_THER',
                           MODELE=MODT,
                           CHAM_MATER=affectMat,
                           )

monSolver = code_aster.PetscSolver( code_aster.Renumbering.Sans )
monSolver.setPreconditioning(code_aster.Preconditioning.Without)

numeDDL = code_aster.ParallelDOFNumbering()
numeDDL.setElementaryMatrix(matr_elem)
numeDDL.computeNumbering()
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")
# numeDDL.debugPrint()

matrAsse = code_aster.AssemblyMatrixTemperatureReal()
matrAsse.appendElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.addDirichletBC(charCine)
matrAsse.build()
test.assertEqual(matrAsse.getType(), "MATR_ASSE_TEMP_R")
# matrAsse.debugPrint()

# retour = vect_elem.assembleVector( numeDDL )
vecass = ASSE_VECTEUR(VECT_ELEM=vect_elem, NUME_DDL=numeDDL)
vcine = CALC_CHAR_CINE(NUME_DDL=numeDDL,
                       CHAR_CINE=charCine,)

# monSolver.matrixFactorization( matrAsse )
# resu = monSolver.solveRealLinearSystemWithDirichletBC(matrAsse, vcine, vecass)
# resu.debugPrint(6)

matrAsse = FACTORISER(reuse=matrAsse,
                      MATR_ASSE=matrAsse,
                      METHODE='PETSC', PRE_COND='JACOBI',
                      )

resu = RESOUDRE(MATR=matrAsse,
                  CHAM_NO=vecass,
                  CHAM_CINE=vcine,
                  ALGORITHME='CR',
                  RESI_RELA=1E-9,)

TEST_RESU(
    CHAM_NO=_F(
        CRITERE='ABSOLU',
        GROUP_NO='Noeud1',
        NOM_CMP='TEMP',
        REFERENCE='AUTRE_ASTER',
        CHAM_GD=resu,
        VALE_CALC=8.33333333333333E-10,
        VALE_REFE=8.33333333333333E-10,
    ))

TEST_RESU(
    CHAM_NO=_F(
        CRITERE='ABSOLU',
        GROUP_NO='Noeud2',
        NOM_CMP='TEMP',
        REFERENCE='AUTRE_ASTER',
        CHAM_GD=resu,
        VALE_CALC=8.33333333333333E-10,
        VALE_REFE=8.33333333333333E-10,
    ))



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

    rank = A.getComm().getRank()
    print('rank=', rank)
    rs, re = A.getOwnershipRange()
    ce, _ = A.getSize()
    rows = N.array(list(range(rs, re)), dtype=petsc4py.PETSc.IntType)
    cols = N.array(list(range(0, ce)), dtype=petsc4py.PETSc.IntType)
    rows = petsc4py.PETSc.IS().createGeneral(rows, comm=A.getComm())
    cols = petsc4py.PETSc.IS().createGeneral(cols, comm=A.getComm())
    (S,) = A.createSubMatrices(rows, cols)
    v = petsc4py.PETSc.Viewer().createASCII("mesh004c_rank"+str(rank)+".out", comm=S.getComm())
    S.view(v)

test.printSummary()

FIN()
