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
import code_aster
from code_aster.Commands import *

code_aster.init("--test")

test = code_aster.TestCase()

rank = code_aster.MPI.ASTER_COMM_WORLD.Get_rank()

MAIL = LIRE_MAILLAGE(FORMAT='MED',
                     PARTITIONNEUR='PTSCOTCH',
                     INFO=1,
                     )
DEFI_GROUP(reuse=MAIL, MAILLAGE=MAIL, CREA_GROUP_NO=_F(TOUT_GROUP_MA='OUI', ))

MATER = DEFI_MATERIAU(ELAS=_F(E=10000.0,
                              NU=0.,
                              RHO=1.0, ),
                      )

affectMat = code_aster.MaterialField(MAIL)
affectMat.addMaterialOnMesh(MATER)
affectMat.build()

MODT = AFFE_MODELE(MAILLAGE=MAIL,
                   AFFE=_F(TOUT='OUI',  # GROUP_MA=('S11',    'S31', 'S12',     'S32'),
                           PHENOMENE='MECANIQUE',
                           MODELISATION='D_PLAN', ), )


CHT1 = AFFE_CHAR_MECA(MODELE=MODT,
                      PRES_REP=_F(TOUT='OUI',
                                  PRES=-1, ), )

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
matr_elem = dComputation.elasticStiffnessMatrix()

monSolver = code_aster.MumpsSolver()

numeDDL = study.getDOFNumbering()
test.assertEqual(numeDDL.getType(), "NUME_DDL_P")

matrAsse = code_aster.AssemblyMatrixDisplacementReal(study)
matrAsse.addElementaryMatrix(matr_elem)
matrAsse.assemble()
test.assertEqual(matrAsse.getType(), "MATR_ASSE_DEPL_R")

matrAsse *= 2.0
matrAsse *= 0.5

matrAsse2 = matrAsse.duplicate()

test.assertNotEqual(matrAsse.getName(), matrAsse2.getName())
test.assertEqual(matrAsse.isEmpty(), matrAsse2.isEmpty())
test.assertEqual(matrAsse.getDOFNumbering(), matrAsse2.getDOFNumbering())
test.assertEqual(matrAsse.getListOfLoads(), matrAsse2.getListOfLoads())
test.assertAlmostEqual(matrAsse.getLagrangeScaling(),
                       matrAsse2.getLagrangeScaling())
petscMat2 = matrAsse2.toPetsc()
ref_matr = 108726.49458749207
test.assertAlmostEqual(ref_matr, petscMat2.norm(), delta=ref_matr*1.e-6)

matrAsse3 = -(-1.0 * matrAsse)
test.assertNotEqual(matrAsse.getName(), matrAsse3.getName())
test.assertEqual(matrAsse.isEmpty(), matrAsse3.isEmpty())
test.assertEqual(matrAsse.getDOFNumbering(), matrAsse3.getDOFNumbering())
test.assertEqual(matrAsse.getListOfLoads(), matrAsse3.getListOfLoads())
test.assertAlmostEqual(matrAsse.getLagrangeScaling(),
                       matrAsse3.getLagrangeScaling())
petscMat3 = matrAsse3.toPetsc()
test.assertAlmostEqual(ref_matr, petscMat3.norm(), delta=ref_matr*1.e-6)

matrAsse4 = 3.*matrAsse + 2.0*matrAsse2 - 4.0*matrAsse3
test.assertNotEqual(matrAsse.getName(), matrAsse4.getName())
test.assertEqual(matrAsse.isEmpty(), matrAsse4.isEmpty())
test.assertEqual(matrAsse.getDOFNumbering(), matrAsse4.getDOFNumbering())
test.assertEqual(matrAsse.getListOfLoads(), matrAsse4.getListOfLoads())
test.assertAlmostEqual(matrAsse.getLagrangeScaling(),
                       matrAsse4.getLagrangeScaling())
petscMat4 = matrAsse4.toPetsc()
test.assertAlmostEqual(ref_matr, petscMat4.norm(), delta=ref_matr*1.e-6)

matrAsse5 = matrAsse.duplicate()
matrAsse5 += matrAsse
matrAsse5 -= matrAsse
test.assertNotEqual(matrAsse.getName(), matrAsse5.getName())
test.assertEqual(matrAsse.isEmpty(), matrAsse5.isEmpty())
test.assertEqual(matrAsse.getDOFNumbering(), matrAsse5.getDOFNumbering())
test.assertEqual(matrAsse.getListOfLoads(), matrAsse5.getListOfLoads())
test.assertAlmostEqual(matrAsse.getLagrangeScaling(),
                       matrAsse5.getLagrangeScaling())
petscMat5 = matrAsse5.toPetsc()
test.assertAlmostEqual(ref_matr, petscMat5.norm(), delta=ref_matr*1.e-6)

monSolver.factorize(matrAsse)

vcine = CALC_CHAR_CINE(NUME_DDL=numeDDL, CHAR_CINE=charCine, )
resu = monSolver.solve(retour, vcine)

ref = 11.55955851173471
test.assertAlmostEqual(ref, resu.norm("NORM_2"), delta=ref*1.e-6)

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
