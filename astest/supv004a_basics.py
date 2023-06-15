# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

import unittest

from code_aster.Utilities import no_new_attributes


class UtilitiesTest(unittest.TestCase):
    """Check for generic utilities."""

    def test01_attrs(self):
        class A:
            a = b = None

            __setattr__ = no_new_attributes(object.__setattr__)

            def __init__(self, a):
                self.a = a
                A.a = a + 2

            @classmethod
            def setB(cls, b):
                cls.b = b

            def bad(self):
                self.c = 0

        o1 = A(1)
        o1.setB(1)
        self.assertEqual(o1.a, 1)
        self.assertEqual(A.a, 3)
        self.assertEqual(o1.b, 1)

        o2 = A(2)
        o2.setB(2)
        self.assertEqual(o2.a, 2)
        self.assertEqual(o2.b, 2)
        self.assertEqual(o1.a, 1)
        self.assertEqual(o1.b, 2)
        self.assertEqual(A.a, 4)

        with self.assertRaises(AttributeError):
            o1.bad()


if __name__ == "__main__":
    unittest.main()
