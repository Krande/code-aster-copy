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

print(repr(user_object))

test.assertIsInstance(user_object.nested, dict)
test.assertIsInstance(user_object.nested["lotod"], list)
test.assertEqual(len(user_object.nested["lotod"]), 2)
test.assertEqual(user_object.nested["lotod"][0], "list")
test.assertIsInstance(user_object.nested["lotod"][1], tuple)
test.assertEqual(len(user_object.nested["lotod"][1]), 2)
test.assertEqual(user_object.nested["lotod"][1][0], "tuple")
test.assertIsInstance(user_object.nested["lotod"][1][1], dict)

mesh = user_object.nested["lotod"][1][1]["dict"]
test.assertIsInstance(mesh, CA.Mesh)
test.assertEqual(mesh.getNumberOfNodes(), 3)

subobject = user_object.subobj
test.assertIsInstance(subobject.attrname, CA.Mesh)
test.assertEqual(subobject.attrname.getNumberOfNodes(), 3)

test.assertEqual(user_object.values[0], "SIEF_ELGA")
test.assertEqual(user_object.values[1], 31)
test.assertAlmostEqual(user_object.values[2], -325.03920740223253)

test.printSummary()

CA.close()
