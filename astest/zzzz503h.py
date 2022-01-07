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

mesh = code_aster.Mesh.buildCube(nrefine=0)

trac = DEFI_FONCTION(
    NOM_PARA="EPSI",
    ABSCISSE=(0.002, 0.1),
    ORDONNEE=(300.0, 300.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

mater = DEFI_MATERIAU(ELAS=_F(E=3.7272000000e10, NU=0.0, RHO=2400.0), TRACTION=_F(SIGM=trac))

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
test.assertEqual(mat.getNumberOfMaterialProperties(), 2, msg="number of material properties")
test.assertEqual(mat.getNumberOfUserMaterialProperties(), 2, msg="number of user material properties")

# test.assertEqual(mat.getNumberOfListOfPropertiesReal(0), 0, msg="list of floats in ELAS")
# test.assertEqual(mat.getNumberOfListOfPropertiesFunction(0), 0, msg="list of functions in ELAS")
# test.assertEqual(mat.getNumberOfListOfPropertiesReal(1), 0, msg="list of floats in TRACTION")
# test.assertEqual(mat.getNumberOfListOfPropertiesFunction(1), 0, msg="list of functions in TRACTION")

# the order of the assignments should be preserved
elas, trac = mat.getVectorOfMaterialProperties()

print("\nchecking ELAS...")
test.assertEqual(elas.getName(), "ELAS", msg="check name")
test.assertEqual(elas.getAsterName(), "ELAS", msg="check name")

test.assertEqual(elas.getNumberOfListOfPropertiesReal(), 0, msg="list of floats in ELAS")
test.assertEqual(elas.getNumberOfListOfPropertiesFunction(), 0, msg="list of functions in ELAS")

test.assertTrue(elas.hasValueReal("E"), msg="Young modulus")
test.assertAlmostEqual(elas.getValueReal("E"), 3.72720e10, msg="Young modulus")
test.assertTrue(elas.hasValueReal("Nu"), msg="Poisson ratio")
test.assertAlmostEqual(elas.getValueReal("Nu"), 0.0, msg="Poisson ratio")
test.assertTrue(elas.hasValueReal("Rho"), msg="Density")
test.assertAlmostEqual(elas.getValueReal("Rho"), 2400.0, msg="Density")

print("\nchecking TRACTION...")
test.assertEqual(trac.getName(), "TRACTION")
test.assertEqual(trac.getAsterName(), "TRACTION")

test.assertEqual(trac.getNumberOfListOfPropertiesReal(), 0, msg="list of floats in TRACTION")
test.assertEqual(trac.getNumberOfListOfPropertiesFunction(), 0, msg="list of functions in TRACTION")

test.assertTrue(trac.hasValueGenericFunction("Sigm"), msg="Traction curve is function")
test.assertFalse(trac.hasValueReal("Sigm"), msg="Traction curve is float")
test.assertFalse(trac.hasValueComplex("Sigm"), msg="Traction curve is complex")
test.assertFalse(trac.hasValueString("Sigm"),msg="Traction curve is string")
test.assertFalse(trac.hasValueTable("Sigm"), msg="Traction curve is table")

test.printSummary()

FIN()
