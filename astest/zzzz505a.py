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


from code_aster.Commands import *
from code_aster import CA
from math import sqrt

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

# For a Mesh
mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

print(mesh.getGroupsOfNodes(), flush=True)

nbNodes = mesh.getNumberOfNodes()

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

dofNume = NUME_DDL(MODELE=model)

elno = CREA_CHAMP(
    TYPE_CHAM="ELGA_DEPL_R",
    NUME_DDL=dofNume,
    OPERATION="AFFE",
    MODELE=model,
    AFFE=_F(TOUT="OUI", NOM_CMP=("DX", "DY", "DZ"), VALE=(-1.0, 2.0, 3.0)),
)

test.assertEqual(elno.getPhysicalQuantity(), "DEPL_R")
test.assertEqual(elno.getFieldType(), "ELGA")
test.assertSequenceEqual(elno.getComponents(), ["DX", "DY", "DZ"])
test.assertEqual(elno.getNumberOfComponents(), 3)

elno_values = elno.getValues()
elno.setValues(elno_values)

TEST_RESU(
    CHAM_ELEM=_F(
        CHAM_GD=elno,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=-1.0,
        VALE_REFE=-1.0,
        NOM_CMP="DX",
        TYPE_TEST="MIN",
    )
)

field = CA.FieldOnNodesReal(mesh, "DEPL_R", {"DX": -1.0, "DY": 2.0, "DZ": 3.0})

TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=field,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=-1.0,
        VALE_REFE=-1.0,
        NOM_CMP="DX",
        TYPE_TEST="MIN",
    )
)

norm_1 = (1.0 + 2.0 + 3.0) * nbNodes
norm_2 = sqrt(14.0 * nbNodes)
norm_inf = 3.0

test.assertEqual(field.size(), 3 * nbNodes)
test.assertSequenceEqual(mesh.getInnerNodes(), range(nbNodes))
test.assertEqual(len(mesh.getInnerNodes()), nbNodes)
test.assertEqual(len(field.getValues()), field.size())
test.assertAlmostEqual(sum([abs(x) for x in field.getValues()]), norm_1)
test.assertAlmostEqual(field.norm("NORM_1"), norm_1)
test.assertAlmostEqual(field.norm("NORM_2"), norm_2)
test.assertAlmostEqual(field.norm("NORM_INFINITY"), norm_inf)
test.assertAlmostEqual(field.norm("NORM_1", ["DX"]), nbNodes)
test.assertAlmostEqual(max(field.getValues()), norm_inf)
test.assertAlmostEqual(field.dot(field), norm_2 * norm_2)

fr = field.restrict(["DX", "DY", "YD"])
test.assertAlmostEqual(fr.norm("NORM_1"), 3.0 * nbNodes)
test.assertAlmostEqual(fr.norm("NORM_1"), fr.toSimpleFieldOnNodes().toFieldOnNodes().norm("NORM_1"))
test.assertEqual(len(fr.getValues(["DY"], ["A"])), 1)
test.assertAlmostEqual(fr.getValues(["DY"], ["A"])[0], 2.0)
test.assertEqual(len(fr.getValues(["DY"], ["ZZ"])), 0)

with test.assertRaises(CA.AsterError):
    ferror = field.restrict(groupsOfNodes=["NOEXI"])

with test.assertRaises(CA.AsterError):
    ferror = field.restrict(["YD"])

f0 = field.copy()
f = field.copy()
f += f
f.getValues()
f2 = 2 * f0
f2.getValues()
f3 = -f + f2
f3.getValues()
test.assertAlmostEqual(f3.norm("NORM_2"), 0)

myField = CA.FieldOnNodesReal(dofNume)
myField.setValues({"DX": 1.0, "DY": 1.0, "DZ": 1.0, "XX": 1234.5})
test.assertAlmostEqual(myField.norm("NORM_1"), 3 * nbNodes)

field2 = field - myField

test.assertAlmostEqual(field2.norm("NORM_2"), sqrt((4.0 + 1.0 + 4.0) * nbNodes))

# For a Parallel Mesh
meshp = LIRE_MAILLAGE(FORMAT="MED", UNITE=20, PARTITIONNEUR="PTSCOTCH")

modelp = AFFE_MODELE(MAILLAGE=meshp, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

nu = NUME_DDL(MODELE=modelp)

fieldp = CREA_CHAMP(
    TYPE_CHAM="NOEU_DEPL_R",
    OPERATION="AFFE",
    MODELE=modelp,
    NUME_DDL=nu,
    AFFE=_F(TOUT="OUI", NOM_CMP=("DX", "DY", "DZ"), VALE=(1.0, 2.0, 3.0)),
)

test.assertEqual(fieldp.size(), 3 * meshp.getNumberOfNodes())
test.assertAlmostEqual(fieldp.norm("NORM_1"), norm_1)
test.assertAlmostEqual(fieldp.norm("NORM_2"), norm_2)
test.assertAlmostEqual(fieldp.norm("NORM_INFINITY"), norm_inf)
test.assertAlmostEqual(fieldp.dot(fieldp), norm_2 * norm_2)

f0 = fieldp.copy()
f = fieldp.copy()
f += f
f.getValuesWithDescription()
f2 = 2 * f0
f2.getValuesWithDescription()
f3 = -f + f2
f3.getValuesWithDescription()
test.assertAlmostEqual(f3.norm("NORM_2"), 0)

# Test TEST_RESU with TEST_TYPE='MIN','MAX', 'SOMME', 'SOMME_ABS'
ftest = fieldp.copy()
mapDOF = ftest.getDescription().getDOFFromNodeAndComponent(False)

values_test = {
    "MAX": 4.0,
    "MIN": -7.0,
    "SOMM": 25 * 1.0 + 4.0 - 7.0,
    "SOMM_ABS": 25 * 1.0 + 4.0 + 7.0,
}
ftest.updateValuePointers()
# for 'MAX' - DX
ftest[mapDOF[1, "DX"]] = 4.0
# for 'MIN' - DX
ftest[mapDOF[15, "DX"]] = -7.0

# Avec NOM_CMP
TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=-7.0,
        VALE_REFE=values_test["MIN"],
        NOM_CMP="DX",
        TYPE_TEST="MIN",
    )
)

TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=4.0,
        VALE_REFE=values_test["MAX"],
        NOM_CMP="DX",
        TYPE_TEST="MAX",
    )
)

TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=22.0,
        VALE_REFE=values_test["SOMM"],
        NOM_CMP="DX",
        TYPE_TEST="SOMM",
    )
)

TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=36.0,
        VALE_REFE=values_test["SOMM_ABS"],
        NOM_CMP="DX",
        TYPE_TEST="SOMM_ABS",
    )
)

# Sans NOM_CMP
TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=-7.0,
        VALE_REFE=values_test["MIN"],
        TYPE_TEST="MIN",
    )
)

TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=4.0,
        VALE_REFE=values_test["MAX"],
        TYPE_TEST="MAX",
    )
)

TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest, REFERENCE="ANALYTIQUE", VALE_CALC=157.0, VALE_REFE=157.0, TYPE_TEST="SOMM"
    )
)

TEST_RESU(
    CHAM_NO=_F(
        CHAM_GD=ftest,
        REFERENCE="ANALYTIQUE",
        VALE_CALC=171.0,
        VALE_REFE=171.0,
        TYPE_TEST="SOMM_ABS",
    )
)

# TEST getValuesWithDescription for FieldOnNodesReal
f_real = CA.FieldOnNodesReal(dofNume)
f_real.setValues(1.0)
vals_real, _ = f_real.getValuesWithDescription("DX")
test.assertEqual(vals_real[0], 1.0)

sf_real = f_real.toSimpleFieldOnNodes()
sf_real_values, sf_real_mask = sf_real.getValues()
test.assertEqual(sf_real_values[0][0], 1.0)
test.assertEqual(sf_real_mask.all(), True)

# TEST getValuesWithDescription for FieldOnNodesComplex
f_complex = CA.FieldOnNodesComplex(dofNume)
f_complex.setValues(1 + 2j)
vals_complex, _ = f_complex.getValuesWithDescription("DX")
test.assertEqual(vals_complex[0], 1 + 2j)

f2_complex = CA.FieldOnNodesComplex(dofNume)
f2_complex.setValues(1 + 5j)
dot_f_f2 = sum(
    f_complex.getValues()[p] * f2_complex.getValues()[p].conjugate()
    for p in range(len(f_complex.getValues(["DX", "DY", "DZ"])))
)
test.assertEqual(f_complex.dot(f2_complex), dot_f_f2)
test.assertEqual(f_complex.norm("NORM_1"), sum(map(abs, f_complex.getValues())))
test.assertEqual(
    f_complex.norm("NORM_2"), sqrt(sum(map(lambda x: abs(x) ** 2, f_complex.getValues())))
)
test.assertEqual(f_complex.norm("NORM_INFINITY"), max(map(abs, f_complex.getValues())))

sf_complex = f_complex.toSimpleFieldOnNodes()
sf_complex_values, sf_complex_mask = sf_complex.getValues()
test.assertEqual(sf_complex_values[0][0], 1 + 2j)
test.assertEqual(sf_complex_mask.all(), True)

#
test.printSummary()
CA.close(INFO_BASE="NON")
