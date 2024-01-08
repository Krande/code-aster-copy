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

CA.init("--test", "--continue", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

# before fixing issue31609, 'pres' was a tuple, 'a' was the list 'pres' and 'b' wad 'a'
test.assertIsInstance(pres, list, type(pres))
test.assertEqual(len(pres), 2, len(pres))
test.assertIsInstance(a, CA.Function, a)
test.assertIsInstance(b, CA.Function, b)

test.assertEqual(pres[0](1.0), 1.0)
test.assertEqual(pres[1][0](1.0), 2.0)
test.assertEqual(a(17), 1.0)
test.assertEqual(b(89), 2.0)


def Y0(X):
    return 0.0


form.setContext(dict(Y0=Y0))
test.assertEqual(form(0.0), 0.0)

test.printSummary()

CA.close()
