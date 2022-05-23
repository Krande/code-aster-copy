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

code_aster.init("--test")

test = code_aster.TestCase()

mesh = code_aster.Mesh.buildCube(refine=0)

trac = DEFI_FONCTION(
    NOM_PARA="EPSI",
    ABSCISSE=(0.002, 0.1),
    ORDONNEE=(300.0, 300.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

cst1 = DEFI_CONSTANTE(VALE=1.0)
cst2 = DEFI_CONSTANTE(VALE=2.0)

mater = DEFI_MATERIAU(
    ELAS=_F(E=3.7272000000e10, NU=0.0, RHO=2400.0),
    TRACTION=_F(SIGM=trac),
    MFRONT=_F(LISTE_COEF=(3.69e10, 0.3, 151.0, 87.0, 2.3)),
    UMAT_FO=_F(LISTE_COEF=(cst1, cst2)),
)

chmat = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=mater))

# check for Material API
print("\nchecking material assignment...")
list_affe = chmat.getVectorOfPartOfMaterialField()
test.assertEqual(len(list_affe), 1, msg="number of zones")

affe = list_affe[0]
list_mate = list_affe[0].getVectorOfMaterial()
test.assertEqual(len(list_mate), 1, msg="number of materials")

mat = list_mate[0]
print("\nchecking material...")
test.assertEqual(mater, mat, msg="same objects")

# MaterialProperty == factor keyword in DEFI_MATERIAU
test.assertEqual(mat.getNumberOfMaterialProperties(), 4, msg="number of material properties")
test.assertEqual(
    mat.getNumberOfUserMaterialProperties(), 2, msg="number of user material properties"
)

# the order of the assignments should be preserved
properties = mat.getVectorOfMaterialProperties()

elas = [occ for occ in properties if occ.getName() == "ELAS"][0]
trac = [occ for occ in properties if occ.getName() == "TRACTION"][0]
mfr = [occ for occ in properties if occ.getName() == "MFRONT"][0]
umat = [occ for occ in properties if occ.getName() == "UMAT_FO"][0]

print("\nchecking ELAS...")
test.assertEqual(elas.getName(), "ELAS", msg="check name")
test.assertEqual(elas.getAsterName(), "ELAS", msg="check name")

test.assertEqual(elas.getNumberOfPropertiesReal(), 4, msg="floats in ELAS")
test.assertEqual(elas.getNumberOfPropertiesComplex(), 0, msg="complex numbers in ELAS")
test.assertEqual(elas.getNumberOfPropertiesString(), 0, msg="strings in ELAS")
test.assertEqual(elas.getNumberOfPropertiesFunction(), 0, msg="functions in ELAS")
test.assertEqual(elas.getNumberOfPropertiesTable(), 0, msg="tables in ELAS")
test.assertEqual(elas.getNumberOfListOfPropertiesReal(), 0, msg="list of floats in ELAS")
test.assertEqual(elas.getNumberOfListOfPropertiesFunction(), 0, msg="list of functions in ELAS")

test.assertTrue(elas.hasValueReal("E"), msg="Young modulus")
test.assertAlmostEqual(elas.getValueReal("E"), 3.72720e10, msg="Young modulus")
test.assertTrue(elas.hasValueReal("Nu"), msg="Poisson ratio")
test.assertAlmostEqual(elas.getValueReal("Nu"), 0.0, msg="Poisson ratio")
test.assertTrue(elas.hasValueReal("Rho"), msg="Density")
test.assertAlmostEqual(elas.getValueReal("Rho"), 2400.0, msg="Density")
test.assertTrue(elas.hasValueReal("Coef_amor"), msg="Damping coefficient")
test.assertAlmostEqual(elas.getValueReal("Coef_amor"), 1.0, msg="Damping coefficient")
test.assertFalse(elas.hasValueGenericFunction("Sigm"), msg="no function")

print("\nchecking TRACTION...")
test.assertEqual(trac.getName(), "TRACTION")
test.assertEqual(trac.getAsterName(), "TRACTION")

test.assertEqual(trac.getNumberOfPropertiesReal(), 0, msg="floats in TRACTION")
test.assertEqual(trac.getNumberOfPropertiesComplex(), 0, msg="complex numbers in TRACTION")
test.assertEqual(trac.getNumberOfPropertiesString(), 0, msg="strings in TRACTION")
test.assertEqual(trac.getNumberOfPropertiesFunction(), 1, msg="functions in TRACTION")
test.assertEqual(trac.getNumberOfPropertiesTable(), 0, msg="tables in TRACTION")
test.assertEqual(trac.getNumberOfListOfPropertiesReal(), 0, msg="list of floats in TRACTION")
test.assertEqual(trac.getNumberOfListOfPropertiesFunction(), 0, msg="list of functions in TRACTION")

test.assertTrue(trac.hasValueGenericFunction("Sigm"), msg="Traction curve is function")
test.assertFalse(trac.hasValueReal("Sigm"), msg="Traction curve is float")
test.assertFalse(trac.hasValueComplex("Sigm"), msg="Traction curve is complex")
test.assertFalse(trac.hasValueString("Sigm"), msg="Traction curve is string")
test.assertFalse(trac.hasValueTable("Sigm"), msg="Traction curve is table")

print("\nchecking MFRONT...")
test.assertEqual(mfr.getName(), "MFRONT")
test.assertEqual(mfr.getAsterName(), "MFRONT")

test.assertEqual(mfr.getNumberOfPropertiesReal(), 0, msg="floats in MFRONT")
test.assertEqual(mfr.getNumberOfPropertiesComplex(), 0, msg="complex numbers in MFRONT")
test.assertEqual(mfr.getNumberOfPropertiesString(), 0, msg="strings in MFRONT")
test.assertEqual(mfr.getNumberOfPropertiesFunction(), 0, msg="functions in MFRONT")
test.assertEqual(mfr.getNumberOfPropertiesTable(), 0, msg="tables in MFRONT")
test.assertEqual(mfr.getNumberOfListOfPropertiesReal(), 1, msg="list of floats in MFRONT")
test.assertEqual(mfr.getNumberOfListOfPropertiesFunction(), 0, msg="list of functions in MFRONT")

print("\nchecking UMAT_FO...")
test.assertEqual(umat.getName(), "UMAT_FO")
test.assertEqual(umat.getAsterName(), "UMAT_FO")

test.assertEqual(umat.getNumberOfPropertiesReal(), 0, msg="floats in UMAT_FO")
test.assertEqual(umat.getNumberOfPropertiesComplex(), 0, msg="complex numbers in UMAT_FO")
test.assertEqual(umat.getNumberOfPropertiesString(), 0, msg="strings in UMAT_FO")
test.assertEqual(umat.getNumberOfPropertiesFunction(), 0, msg="functions in UMAT_FO")
test.assertEqual(umat.getNumberOfPropertiesTable(), 0, msg="tables in UMAT_FO")
test.assertEqual(umat.getNumberOfListOfPropertiesReal(), 0, msg="list of floats in UMAT_FO")
test.assertEqual(umat.getNumberOfListOfPropertiesFunction(), 1, msg="list of functions in UMAT_FO")

test.printSummary()

FIN()
