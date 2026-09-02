# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
from code_aster.visu.Pairing.convexPointSet import (
    ConvexPointSet,
    CSStrat2DCodim1,
    CSStrat2DCodim0,
    CSStrat1D,
)
import numpy as np


class TestConvexSetPoint(unittest.TestCase):
    """Tests the ConvexSetPoint class"""

    def test_init_ConvexSetPoint(self):
        set_1d = ConvexPointSet(dim=1, points=[[0], [1]], codim=0)
        self.assertIsInstance(set_1d._strategy, CSStrat1D)

        set_2d0 = ConvexPointSet(dim=2, points=[[0, 0], [1, 0], [0, 1]], codim=0)
        self.assertIsInstance(set_2d0._strategy, CSStrat2DCodim0)

        set_2d1 = ConvexPointSet(dim=2, points=[[0, 0, 0], [1, 0, 0], [0, 1, 0]], codim=1)
        self.assertIsInstance(set_2d1._strategy, CSStrat2DCodim1)


class TestCSStrat1D(unittest.TestCase):
    """Unit tests for the 1D strategy (CSStrat1D)."""

    def test_area_is_length(self):
        """Test that the 'area' of a 1D set is its length."""
        points = np.array([[5.0, 2.0], [1.0, 2.0]])
        strategy = CSStrat1D(points)
        self.assertAlmostEqual(strategy.computeArea(), 4.0)

    def test_barycenter(self):
        """Test the barycenter calculation for a 1D set."""
        points = np.array([[5.0, 2.0], [1.0, 2.0]])
        strategy = CSStrat1D(points)
        expected_bary = np.array([3.0, 2.0])
        np.testing.assert_allclose(strategy.computeBary(), expected_bary)


class TestCSStrat2DCodim0(unittest.TestCase):
    """Unit tests for the 2D, Codim 0 strategy (CSStrat2DCodim0)."""

    def setUp(self):
        # - A 2x2 square, with points shuffled to test sorting
        self.shuffled_points = np.array([[0, 0], [2, 2], [2, 0], [0, 2]])

    def test_area_of_square(self):
        """Test the area of a 2x2 square in 2D."""
        strategy = CSStrat2DCodim0(self.shuffled_points)
        self.assertAlmostEqual(strategy.computeArea(), 4.0)

    def test_barycenter_of_square(self):
        """Test the barycenter of a 2x2 square in 2D."""
        strategy = CSStrat2DCodim0(self.shuffled_points)
        expected_bary = np.array([1.0, 1.0])
        np.testing.assert_allclose(strategy.computeBary(), expected_bary)


class TestCSStrat2DCodim1(unittest.TestCase):
    """Unit tests for the 2D, Codim 1 strategy."""

    def setUp(self):
        # - A 2x2 square on the XY plane, with points shuffled
        self.shuffled_points_xy = np.array([[0, 0, 0], [2, 2, 0], [2, 0, 0], [0, 2, 0]])
        # - A 2x2 square on the XZ plane (tilted)
        self.tilted_points_xz = np.array([[0, 0, 0], [2, 0, 0], [2, 0, 2], [0, 0, 2]])

    def test_area_of_square_in_3d(self):
        """Test the area of a 2x2 square embedded in 3D."""
        strategy = CSStrat2DCodim1(self.shuffled_points_xy)
        self.assertAlmostEqual(strategy.computeArea(), 4.0)

    def test_barycenter_of_square_in_3d(self):
        """Test the barycenter of a 2x2 square embedded in 3D."""
        strategy = CSStrat2DCodim1(self.shuffled_points_xy)
        expected_bary = np.array([1.0, 1.0, 0.0])
        np.testing.assert_allclose(strategy.computeBary(), expected_bary)

    def test_area_of_tilted_square(self):
        """Test area calculation for a polygon not aligned with major planes."""
        strategy = CSStrat2DCodim1(self.tilted_points_xz)
        self.assertAlmostEqual(strategy.computeArea(), 4.0)

    def test_local_basis_computation(self):
        """Test that the local basis is orthonormal."""
        strategy = CSStrat2DCodim1(self.tilted_points_xz)
        basis = strategy.computeLocalBasis()
        vec1, vec2, norm = basis[0], basis[1], basis[2]
        # - Check for orthogonality
        self.assertAlmostEqual(np.dot(vec1, vec2), 0)
        self.assertAlmostEqual(np.dot(vec1, norm), 0)
        # - Check for unit length
        self.assertAlmostEqual(np.linalg.norm(vec1), 1)
        self.assertAlmostEqual(np.linalg.norm(vec2), 1)
        self.assertAlmostEqual(np.linalg.norm(norm), 1)


if __name__ == "__main__":
    unittest.main()
