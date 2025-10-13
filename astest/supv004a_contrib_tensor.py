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

import unittest

import numpy as NP

from code_aster.MacroCommands.Contrib.TensorModule import (
    Tensor,
    div,
    flatten,
    grad,
    gradsym,
    isTensor,
    laplacien,
)

from code_aster.MacroCommands.Contrib.HookeTensor import HookeOrthotropic

try:
    from code_aster.Utilities import sympy

    X, Y, Z = sympy.symbols("X Y Z")
    ASTER_HAVE_SYMPY = True
except ImportError:
    ASTER_HAVE_SYMPY = False


class TensorUnitTest(unittest.TestCase):
    def setUp(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.U = Tensor(
            NP.array(
                (
                    [X**3, sympy.sin(X), sympy.exp(X)],
                    [Y**3, sympy.sin(Y), sympy.exp(Y)],
                    [Z**3, sympy.sin(Z), sympy.exp(Z)],
                )
            )
        )

    def testType(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(isTensor(self.U), 1)

    def testRank(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(self.U.rank, 2)
        self.assertEqual(grad(self.U).rank, 3)

    def testGrad(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(
            grad(self.U),
            Tensor(
                NP.array(
                    [
                        [[3 * X**2, 0, 0], [0, 3 * Y**2, 0], [0, 0, 3 * Z**2]],
                        [[sympy.cos(X), 0, 0], [0, sympy.cos(Y), 0], [0, 0, sympy.cos(Z)]],
                        [[sympy.exp(X), 0, 0], [0, sympy.exp(Y), 0], [0, 0, sympy.exp(Z)]],
                    ]
                )
            ),
        )

    def testGradSym(self):
        if not ASTER_HAVE_SYMPY:
            return
        # attention: sensible 3.0*X**2 != 3*X**2
        self.assertEqual(
            gradsym(self.U),
            Tensor(
                NP.array(
                    [
                        [
                            [3.0 * X**2, 0.5 * sympy.cos(X), 0.5 * sympy.exp(X)],
                            [0, 1.5 * Y**2, 0],
                            [0, 0, 1.5 * Z**2],
                        ],
                        [
                            [0.5 * sympy.cos(X), 0, 0],
                            [1.5 * Y**2, 1.0 * sympy.cos(Y), 0.5 * sympy.exp(Y)],
                            [0, 0, 0.5 * sympy.cos(Z)],
                        ],
                        [
                            [0.5 * sympy.exp(X), 0, 0],
                            [0, 0.5 * sympy.exp(Y), 0],
                            [1.5 * Z**2, 0.5 * sympy.cos(Z), 1.0 * sympy.exp(Z)],
                        ],
                    ]
                )
            ),
        )

    def testLaplacien(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(
            laplacien(self.U),
            Tensor(
                NP.array(
                    [
                        [6 * X, 6 * Y, 6 * Z],
                        [-sympy.sin(X), -sympy.sin(Y), -sympy.sin(Z)],
                        [sympy.exp(X), sympy.exp(Y), sympy.exp(Z)],
                    ]
                )
            ),
        )

    def testDivergence(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(
            div(grad(self.U)),
            Tensor(
                NP.array(
                    [
                        [6 * X, 6 * Y, 6 * Z],
                        [-sympy.sin(X), -sympy.sin(Y), -sympy.sin(Z)],
                        [sympy.exp(X), sympy.exp(Y), sympy.exp(Z)],
                    ]
                )
            ),
        )

    def testProduitDoubleContracte(self):
        if not ASTER_HAVE_SYMPY:
            return
        TensO4Sym = Tensor(
            NP.array(
                [
                    [
                        [[400.0, 0.0, 0.0], [0.0, 200.0, 0.0], [0.0, 0.0, 200.0]],
                        [[0.0, 66.66666667, 0.0], [66.66666667, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 133.33333333], [0.0, 0.0, 0.0], [133.33333333, 0.0, 0.0]],
                    ],
                    [
                        [[0.0, 66.66666667, 0.0], [66.66666667, 0.0, 0.0], [0.0, 0.0, 0.0]],
                        [[200.0, 0.0, 0.0], [0.0, 233.33333333, 0.0], [0.0, 0.0, 166.66666667]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 66.66666667], [0.0, 66.66666667, 0.0]],
                    ],
                    [
                        [[0.0, 0.0, 133.33333333], [0.0, 0.0, 0.0], [133.33333333, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0], [0.0, 0.0, 66.66666667], [0.0, 66.66666667, 0.0]],
                        [[200.0, 0.0, 0.0], [0.0, 166.66666667, 0.0], [0.0, 0.0, 233.33333333]],
                    ],
                ]
            )
        )

        self.assertEqual(
            TensO4Sym.produitDoubleContracte(self.U),
            Tensor(
                NP.array(
                    [
                        [
                            200.000000000000 * sympy.sin(Y)
                            + 400.000000000000 * X**3
                            + 200.000000000000 * sympy.exp(Z),
                            66.6666666700000 * sympy.sin(X) + 66.6666666700000 * Y**3,
                            133.333333330000 * Z**3 + 133.333333330000 * sympy.exp(X),
                        ],
                        [
                            66.6666666700000 * sympy.sin(X) + 66.6666666700000 * Y**3,
                            233.333333330000 * sympy.sin(Y)
                            + 200.000000000000 * X**3
                            + 166.666666670000 * sympy.exp(Z),
                            66.6666666700000 * sympy.sin(Z) + 66.6666666700000 * sympy.exp(Y),
                        ],
                        [
                            133.333333330000 * Z**3 + 133.333333330000 * sympy.exp(X),
                            66.6666666700000 * sympy.sin(Z) + 66.6666666700000 * sympy.exp(Y),
                            166.666666670000 * sympy.sin(Y)
                            + 200.000000000000 * X**3
                            + 233.333333330000 * sympy.exp(Z),
                        ],
                    ]
                )
            ),
        )

    def testproduitSimpleContracte(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(
            self.U.produitSimpleContracte(Tensor(NP.array([-1, 0, 0]))),
            Tensor(NP.array([-(X**3), -(Y**3), -(Z**3)])),
        )


class TensorUnitTest2(unittest.TestCase):
    def setUp(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.U = Tensor(NP.array(([X**3, Y**3, Z**3])))

    def testType(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(isTensor(self.U), 1)

    def testRank(self):
        if not ASTER_HAVE_SYMPY:
            return
        self.assertEqual(self.U.rank, 1)
        self.assertEqual(grad(self.U).rank, 2)

    def testProduitDoubleContracte(self):
        if not ASTER_HAVE_SYMPY:
            return
        tensDiff = HookeOrthotropic(
            200.0, 100.0, 150.0, 0.4, 0.2, 0.3, 100.0, 100.0, 200.0
        ).produitDoubleContracte(Tensor(NP.ones((3, 3)))) - Tensor(
            NP.array(
                [
                    [375.52155772, 200.0, 200.0],
                    [200.0, 273.99165508, 400.0],
                    [200.0, 400.0, 329.62447844],
                ]
            )
        )
        diff = max(NP.fabs(flatten(tensDiff.array.tolist())))
        self.assertAlmostEqual(diff, 0.0, 8)


if __name__ == "__main__":
    unittest.main()
