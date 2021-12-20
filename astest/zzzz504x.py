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

rank = code_aster.MPI.COMM_WORLD.Get_rank()

MAIL = LIRE_MAILLAGE(FORMAT='MED',
                     PARTITIONNEUR='PTSCOTCH',
                     INFO=1,
                     )
DEFI_GROUP(reuse=MAIL, MAILLAGE=MAIL, CREA_GROUP_NO=_F(TOUT_GROUP_MA='OUI', ))

MATER = DEFI_MATERIAU(ELAS=_F(E=10000.0,
                              NU=0.,
                              RHO=1.0, ),
                      );

affectMat = code_aster.MaterialField(MAIL)
affectMat.addMaterialsOnMesh(MATER)
affectMat.buildWithoutExternalStateVariables()

MODT = AFFE_MODELE(MAILLAGE=MAIL,
                   AFFE=_F(TOUT='OUI',  # GROUP_MA=('S11',    'S31', 'S12',     'S32'),
                           PHENOMENE='MECANIQUE',
                           MODELISATION='D_PLAN', ), )

# Ne fonctionne à voir pourquoi ?
# CHT1 = AFFE_CHAR_MECA(MODELE=MODT,
#                      PESANTEUR=_F(GROUP_MA='Bas1', GRAVITE=1.0,
#                                   DIRECTION=(0.0, -1.0, 0.0),),
#                      INFO=1,
#                      VERI_NORM='NON',)

CHT1 = AFFE_CHAR_MECA(MODELE=MODT,
                      PRES_REP=_F(TOUT='OUI',
                                  PRES=-1, ), );

charCine = code_aster.MechanicalDirichletBC(MODT)
charCine.addBCOnCells(code_aster.PhysicalQuantityComponent.Dx, 0., "Bas1")
charCine.addBCOnCells(code_aster.PhysicalQuantityComponent.Dy, 1., "Bas1")
charCine.addBCOnCells(code_aster.PhysicalQuantityComponent.Dx, 0., "Bas3")
charCine.addBCOnCells(code_aster.PhysicalQuantityComponent.Dy, -1., "Bas3")
charCine.build()

study = code_aster.PhysicalProblem(MODT, affectMat)
study.addDirichletBC(charCine)
study.addLoad(CHT1)
study.computeDOFNumbering()
dComputation = code_aster.DiscreteComputation(study)
# compute Neumann
retour = dComputation.neumann([0, 0, 0])
matr_elem = dComputation.computeMechanicalStiffnessMatrix()

monSolver = code_aster.MumpsSolver()

numeDDL = code_aster.ParallelDOFNumbering()
numeDDL.setElementaryMatrix(matr_elem)
numeDDL.computeNumbering()
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")

matrAsse = code_aster.AssemblyMatrixDisplacementReal()
matrAsse.appendElementaryMatrix(matr_elem)
matrAsse.setDOFNumbering(numeDDL)
matrAsse.addDirichletBC(charCine)
matrAsse.build()
test.assertEqual(matrAsse.getType(), "MATR_ASSE_DEPL_R")


monSolver.factorize(matrAsse)

vcine = CALC_CHAR_CINE(NUME_DDL=numeDDL, CHAR_CINE=charCine, )
resu = monSolver.solveWithDirichletBC(retour, vcine)

TEST_RESU(
    CHAM_NO=_F(
        CRITERE='ABSOLU',
        GROUP_NO='Point1',
        NOM_CMP='DY',
        REFERENCE='AUTRE_ASTER',
        CHAM_GD=resu,
        VALE_CALC=1.2743615988908403,
        VALE_REFE=1.2743615988908403,
    ))

TEST_RESU(
    CHAM_NO=_F(
        CRITERE='ABSOLU',
        GROUP_NO='Point4',
        NOM_CMP='DY',
        REFERENCE='AUTRE_ASTER',
        CHAM_GD=resu,
        VALE_CALC=-1.274361598890839,
        VALE_REFE=-1.274361598890839,
    ))

test.printSummary()

FIN()
