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
from math import sqrt

code_aster.init("--test")

test = code_aster.TestCase()

# For a Mesh
mesh = LIRE_MAILLAGE(FORMAT='MED', UNITE=20)

nbNodes = mesh.getNumberOfNodes()

model = AFFE_MODELE(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                            PHENOMENE='MECANIQUE',
                            MODELISATION='3D', ), )

dofNume = NUME_DDL(MODELE=model, )

field = CREA_CHAMP(TYPE_CHAM='NOEU_DEPL_R',
                   OPERATION='AFFE',
                   NUME_DDL=dofNume,
                   MODELE=model,
                   AFFE=_F(TOUT='OUI',
                           NOM_CMP=('DX', 'DY', 'DZ'),
                           VALE=(-1.0, 2.0, 3.0),
                           ), )

TEST_RESU(CHAM_NO=_F(CHAM_GD=field,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=-1.0,
                     VALE_REFE=-1.0,
                     NOM_CMP='DX',
                     TYPE_TEST='MIN', )
          )

norm_1 = (1.0 + 2.0 + 3.0) * nbNodes
norm_2 = sqrt(14.0 * nbNodes)
norm_inf = 3.0

test.assertEqual(field.size(), 3 * nbNodes)
test.assertSequenceEqual(mesh.getInnerNodes(), range(1, nbNodes + 1))
test.assertEqual(len(mesh.getInnerNodes()), nbNodes)
test.assertEqual(len(field.getValues()), field.size())
test.assertAlmostEqual(sum([abs(x) for x in field.getValues()]), norm_1)
test.assertAlmostEqual(field.norm("NORM_1"), norm_1)
test.assertAlmostEqual(field.norm("NORM_2"), norm_2)
test.assertAlmostEqual(field.norm("NORM_INFINITY"), norm_inf)
test.assertAlmostEqual(max(field.getValues()), norm_inf)
test.assertAlmostEqual(field.dot(field), norm_2 * norm_2)

f0 = field.duplicate()
f = field.duplicate()
f += f
f.EXTR_COMP().valeurs
f2 = 2 * f0
f2.EXTR_COMP().valeurs
f3 = -f + f2
f3.EXTR_COMP().valeurs
test.assertAlmostEqual(f3.norm("NORM_2"), 0)

myField = code_aster.FieldOnNodesReal(dofNume)
myField.setValues(1.0)
test.assertAlmostEqual(myField.norm("NORM_1"), 3 * nbNodes)

field2 = field - myField

test.assertAlmostEqual(field2.norm("NORM_2"), sqrt((4.0 + 1.0 + 4.0) * nbNodes))

# For a Parallel Mesh
meshp = LIRE_MAILLAGE(FORMAT='MED', UNITE=20, PARTITIONNEUR="PTSCOTCH")

modelp = AFFE_MODELE(MAILLAGE=meshp,
                     AFFE=_F(TOUT='OUI',
                             PHENOMENE='MECANIQUE',
                             MODELISATION='3D', ), )

fieldp = CREA_CHAMP(TYPE_CHAM='NOEU_DEPL_R',
                    OPERATION='AFFE',
                    MODELE=modelp,
                    AFFE=_F(TOUT='OUI',
                            NOM_CMP=('DX', 'DY', 'DZ'),
                            VALE=(1.0, 2.0, 3.0),
                            ), )

test.assertEqual(fieldp.size(), 3 * meshp.getNumberOfNodes())
test.assertAlmostEqual(fieldp.norm("NORM_1"), norm_1)
test.assertAlmostEqual(fieldp.norm("NORM_2"), norm_2)
test.assertAlmostEqual(fieldp.norm("NORM_INFINITY"), norm_inf)
test.assertAlmostEqual(fieldp.dot(fieldp), norm_2 * norm_2)

f0 = fieldp.duplicate()
f = fieldp.duplicate()
f += f
f.EXTR_COMP().valeurs
f2 = 2 * f0
f2.EXTR_COMP().valeurs
f3 = -f + f2
f3.EXTR_COMP().valeurs
test.assertAlmostEqual(f3.norm("NORM_2"), 0)

# Test TEST_RESU with TEST_TYPE='MIN','MAX', 'SOMME', 'SOMME_ABS'

test.printSummary()
ftest = fieldp.duplicate()
values_test = {'MAX': 4.0, 'MIN': -7.0, 'SOMM': 25 * 1.0 + 4.0 - 7.0, 'SOMM_ABS': 25 * 1.0 + 4.0 + 7.0}
# for 'MAX' - DX
ftest[0] = 4.0
# for 'MIN' - DX
ftest[3] = -7.0

# Avec NOM_CMP
TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=-7.0,
                     VALE_REFE=values_test['MIN'],
                     NOM_CMP='DX',
                     TYPE_TEST='MIN', )
          )

TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=4.0,
                     VALE_REFE=values_test['MAX'],
                     NOM_CMP='DX',
                     TYPE_TEST='MAX', )
          )

TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=22.0,
                     VALE_REFE=values_test['SOMM'],
                     NOM_CMP='DX',
                     TYPE_TEST='SOMM', )
          )

TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=36.0,
                     VALE_REFE=values_test['SOMM_ABS'],
                     NOM_CMP='DX',
                     TYPE_TEST='SOMM_ABS', )
          )

# Sans NOM_CMP
TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=-7.0,
                     VALE_REFE=values_test['MIN'],
                     TYPE_TEST='MIN', )
          )

TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=4.0,
                     VALE_REFE=values_test['MAX'],
                     TYPE_TEST='MAX', )
          )

TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=157.0,
                     VALE_REFE=157.0,
                     TYPE_TEST='SOMM', )
          )

TEST_RESU(CHAM_NO=_F(CHAM_GD=ftest,
                     REFERENCE='ANALYTIQUE',
                     VALE_CALC=171.0,
                     VALE_REFE=171.0,
                     TYPE_TEST='SOMM_ABS', )
          )

# TEST EXTR_COMP for FieldOnNodesComplex

f_complex = code_aster.FieldOnNodesComplex(dofNume)
f_complex.setValues(1+2j)
vals_complex = f_complex.EXTR_COMP('DX').valeurs
test.assertEqual(vals_complex[0], 1+2j)

code_aster.close()
