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

CA.init("--test", "--continue")

test = CA.TestCase()

test.assertEqual(matcomp0.size(), 2, msg="number of material properties")
test.assertCountEqual(matcomp0.getMaterialNames(), ["ELAS", "VISCOCHAB"])

f_alpha = matcomp0.getFunction("ELAS", "ALPHA")
test.assertEqual(f_alpha.userName, "ALPH")
test.assertEqual(f_alpha, ALPH)

test.assertEqual(matcmplx.size(), 1, msg="number of material properties")
test.assertCountEqual(matcmplx.getMaterialNames(), ["ELAS_VISCO"])

test.assertEqual(matcmplx.getValueReal("ELAS_VISCO", "RHO"), 1460.0, msg="get value of RHO")
test.assertEqual(matcmplx.getValueComplex("ELAS_VISCO", "G"), G, msg="get value of G")
test.assertEqual(matcmplx.getValueComplex("ELAS_VISCO", "NU"), nu, msg="get value of NU")

with test.assertRaisesRegex(RuntimeError, "property not found"):
    matcmplx.getValueReal("ELAS_VISCO", "G")

# check for copy constructor
copied = CA.Material(matthm)

test.assertEqual(copied.size(), 7, msg="number of material properties")
test.assertCountEqual(
    copied.getMaterialNames(),
    ["ELAS", "CJS", "THM_LIQU", "THM_GAZ", "THM_VAPE_GAZ", "THM_DIFFU", "THM_INIT"],
)
fLIQ = copied.getFunction("THM_LIQU", "VISC")
test.assertEqual(fLIQ.userName, "VISCOLIQ")
test.assertEqual(copied.getValueReal("THM_LIQU", "RHO"), 1000.0, msg="get value of RHO")
with test.assertRaisesRegex(RuntimeError, "property not found"):
    copied.getValueComplex("THM_LIQU", "BETA")

test.printSummary()

CA.close()
