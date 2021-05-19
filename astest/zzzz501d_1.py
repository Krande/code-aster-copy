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

from zzzz501d_user import ComplexUserObject

code_aster.init("--test", "--continue")
test = code_aster.TestCase()

values = ["SIEF_ELGA", 31, -325.03920740223253]

user_object = ComplexUserObject(MA, U2, values)
print(repr(user_object))

test.assertEqual(user_object.values[0], "SIEF_ELGA")
test.assertEqual(user_object.values[1], 31)
test.assertAlmostEqual(user_object.values[2], -325.03920740223253)

test.printSummary()

code_aster.close()
