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

from math import pi
import numpy as np


from code_aster.Commands import *
from code_aster import CA

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

fsin = CA.Function()
fsin.setParameterName("INST")
fsin.setResultName("TEMP")
fsin.setInterpolation("LIN LIN")
test.assertEqual(fsin.getType(), "FONCTION_SDASTER", msg="fsin: check type")

with test.assertRaises(RuntimeError, msg="fsin: interp type"):
    fsin.setInterpolation("invalid")

fsin.setExtrapolation("CC")

# check properties assignment
prop = fsin.getProperties()
test.assertEqual(prop[1:5], ["LIN LIN", "INST", "TEMP", "CC"], msg="fsin: properties")

# values assignment
n = 12
valx = np.arange(n) * 2.0 * pi / n
valy = np.zeros(n)

# sizes checking
fsin.setValues(valx, valy)
with test.assertRaisesRegex(RuntimeError, "length.*be equal"):
    fsin.setValues(valx, [0.0, 1.0])
with test.assertRaisesRegex(RuntimeError, "function size is"):
    fsin.setValues([0.0, 1.0], [0.0, 1.0])

# assign correct values
valy = np.sin(valx) * -1.0
fsin.setValues(valx, valy)

# fsin.debugPrint(6)

test.assertAlmostEqual(fsin(pi / 2.0), -1.0, msg="fsin(pi/2)")

# check Function.abs()
fabs = fsin.abs()
arrabs = fabs.getValuesAsArray()
test.assertTrue(np.all(arrabs[:, 1] >= 0.0), msg="fsin: abs values")
test.assertAlmostEqual(fabs(pi / 2.0), 1.0, msg="fabs(pi/2)")

values = fsin.getValuesAsArray()
test.assertEqual(values.shape, (n, 2), msg="fsin: shape")

# complex
fcmpl = CA.FunctionComplex()
fcmpl.setParameterName("INST")
fcmpl.setResultName("TEMP")
fcmpl.setInterpolation("LIN LIN")
test.assertEqual(fcmpl.getType(), "FONCTION_C", msg="fcmpl: type")

valz = np.zeros(2 * n)
fcmpl.setValues(valx, valz)

with test.assertRaisesRegex(RuntimeError, "length.*ordinates.*twice.*absc"):
    fcmpl.setValues(valx, [0.0, 1.0])
with test.assertRaisesRegex(RuntimeError, "function size is"):
    fcmpl.setValues([0.0, 1.0], [0.0, 1.0, 0.1, 1.2])

valz = np.vstack([valy, valy]).transpose().ravel()
fcmpl.setValues(valx, valz)

test.assertAlmostEqual(fcmpl(pi / 2.0), -1.0 - 1.0j, msg="fcmpl(pi/2)")

DF1 = DEFI_FONCTION(
    NOM_PARA="INST",
    NOM_RESU="DEPL",
    VERIF="CROISSANT",
    PROL_DROITE="LINEAIRE",
    ABSCISSE=(0.0, 1.0, 2.0),
    ORDONNEE=(0.0, 1.0, 3.0),
)

DF2 = DEFI_FONCTION(
    NOM_PARA="INST",
    NOM_RESU="DEPL",
    INTERPOL="LOG",
    PROL_GAUCHE="LINEAIRE",
    VALE=(3.0, 3.0, 4.0, 4.0, 5.0, 5.0),
)
# DF1.debugPrint()

DN1 = DEFI_NAPPE(
    NOM_PARA="AMOR",
    NOM_RESU="ACCE",
    VERIF="CROISSANT",
    INTERPOL="LOG",
    PROL_DROITE="CONSTANT",
    PROL_GAUCHE="CONSTANT",
    PARA=(0.01, 0.02),
    FONCTION=(DF1, DF2),
)

values2 = DN1.getValues()
nbv2 = len(values2[1]) // 2
test.assertEqual(values2[1][nbv2:], [3.0, 4.0, 5.0], msg="DN1: values")

parameters = DN1.getParameters()
test.assertEqual(parameters, [0.01, 0.02], msg="DN1: parameters")

# DN1.debugPrint()
test.printSummary()

FIN()
