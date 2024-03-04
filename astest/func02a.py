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

test = CA.TestCase()

CA.init("--test", debug=True, ERREUR=_F(ALARME="EXCEPTION"))

one = 1

two = DEFI_CONSTANTE(VALE=2.0)

form_two = FORMULE(NOM_PARA="INST", VALE="two(INST)", two=two)


def SIYYp(x, y, z):
    a0 = 2.0
    ax = -1.0
    ay = -3.0
    az = 1.0
    return (a0 + ax * x + ay * y + az * z) * 1.0e6


SIYY = FORMULE(VALE="SIYYp(X,Y,Z)", SIYYp=SIYYp, NOM_PARA=["X", "Y", "Z"])

test.assertEqual(SIYY(0, 0, 0), 2.0e6)
test.assertEqual(two(0), 2.0)
test.assertEqual(form_two(0), 2.0)

CA.saveObjects()

test.printSummary()

FIN()
