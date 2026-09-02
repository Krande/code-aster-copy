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

from code_aster.visu.Pairing.pairingObjects import (
    PairingObject,
    PairingAnalysisAsterFromPkl,
    PairingAnalysisAster,
)
import numpy as np


## - TEST FOR CLASSES IN convexPointSet
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


## - TEST FOR CLASSES IN pairingObjects


class TestPairingObject(unittest.TestCase):

    def setUp(self):
        self.coords_test = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1]])
        self.connect_test = [[0, 1], [1, 2], [2, 3], [1, 4], [4, 5], [5, 2]]
        self.list_pairs_test = np.array([[0, 1], [0, 2], [3, 4]])

    class DummyPairingObject(PairingObject):
        def setMeshInfos(self, *arg, **kwargs):
            self._flag_MeshInfos = True

        def setPairingInfos(self, *arg, **kwargs):
            self._flag_PairingInfos = True

        def setCellInfos(self, *arg, **kwargs):
            self._flag_CellsInfos = True

    def test_initError(self):
        """Test the abstract instanciation of PairingObject"""
        with self.assertRaises(TypeError):
            PairingObject(2, "m", "mi", "s", "si")

    def test_instanciationMeshInformationError(self):
        """Test the getNodesCoordFromCellIndices logic"""
        obj = self.DummyPairingObject(2, "m", "mi", "s", "si")
        # Flag _flag_MeshInfos should be False by default
        with self.assertRaisesRegex(ValueError, "Mesh informations have not been implemented."):
            obj.getNodesCoordsFromCellIndices([0, 1])

    def test_instanciationMeshInformation(self):
        """Test the getNodesCoordFromCellIndices logic"""
        obj = self.DummyPairingObject(2, "m", "mi", "s", "si")

        obj._coords = self.coords_test
        obj._asterConnectivity = self.connect_test
        obj._flag_MeshInfos = True

        cell_indices = [0, 3]
        expected_node_indices = [0, 1, 4]
        expected_coords = self.coords_test[expected_node_indices]

        coords, node_indices = obj.getNodesCoordsFromCellIndices(cell_indices)

        self.assertListEqual(node_indices, expected_node_indices)
        np.testing.assert_array_equal(coords, expected_coords)

    def test_computebasicInfosFromPairs_logic(self):
        """Test the computebasicInfosFromPairs logic"""
        obj = self.DummyPairingObject(2, "m", "mi", "s", "si")
        obj._listPairs = self.list_pairs_test

        obj.computebasicInfosFromPairs()

        expected_basic_info = np.array([[0, 2], [3, 1]])
        np.testing.assert_array_equal(obj._listPairsBasicInfo, expected_basic_info)

        expected_dict = {
            0: {"indicesCell": [1, 2], "indicesPairs": [0, 1]},
            3: {"indicesCell": [4], "indicesPairs": [2]},
        }
        self.assertDictEqual(obj._listPairsDict, expected_dict)

    def test_getSlaveCellsPaired_logic(self):
        """Test the getSlaveCellsPaired logic"""
        obj = self.DummyPairingObject(2, "m", "mi", "s", "si")
        obj._listPairs = self.list_pairs_test

        self.assertIsNone(obj._listPairsBasicInfo)

        slave_cells = obj.getSlaveCellsPaired()

        self.assertIsNotNone(obj._listPairsBasicInfo)
        np.testing.assert_array_equal(slave_cells, np.array([0, 3]))


class TestPairingAnalysisAsterFromPkl(unittest.TestCase):
    def test_init_PairingAnalysisAsterFromPkl(self):
        """Verifies that the constructor correctly initializes the object."""
        # - Prepare the constructor
        coords = np.array([[0.0, 0.0], [1.0, 1.0]])
        connectivity = [[0, 1]]
        listPairs = np.array([[0, 1]])
        listIntersectionPts = np.array([[0.5, 0.5]])
        listQuadraturePts = np.array([[0.25, 0.25]])
        indices_slv = [0, 2, 4]
        indices_mas = [1, 3, 5]
        indices_do_slv = [10, 12, 14]
        indices_do_mas = [11, 13, 15]

        obj = PairingAnalysisAsterFromPkl(
            dimension=2,
            masterDomain="m",
            masterInterface="mi",
            slaveDomain="s",
            slaveInterface="si",
        )

        self.assertFalse(obj._flag_MeshInfos, "_flag_MeshInfos should be False")
        self.assertFalse(obj._flag_PairingInfos, "_flag_PairingInfos should be False")
        self.assertFalse(obj._flag_CellsInfos, "_flag_CellsInfos should be False")

        obj.setMeshInfos(coords, connectivity)
        self.assertTrue(obj._flag_MeshInfos, "_flag_MeshInfos should be True")

        # - Check mesh information (from setMeshInfos)
        np.testing.assert_array_equal(
            obj._coords, coords, "_coords attribute was not initialized correctly."
        )
        self.assertEqual(
            obj._asterConnectivity,
            connectivity,
            "_asterConnectivity attribute was not initialized correctly.",
        )

        obj.setPairingInfos(listPairs, listIntersectionPts, listQuadraturePts)
        self.assertTrue(obj._flag_MeshInfos, "_flag_MeshInfos should be True")

        # - Check pairing information (from setPairingInfos)
        np.testing.assert_array_equal(
            obj._listPairs, listPairs, "_listPairs attribute was not initialized correctly."
        )
        np.testing.assert_array_equal(
            obj._listIntersectionPts,
            listIntersectionPts,
            "_listIntersectionPts attribute was not initialized correctly.",
        )
        np.testing.assert_array_equal(
            obj._listQuadraturePts,
            listQuadraturePts,
            "_listQuadraturePts attribute was not initialized correctly.",
        )

        obj.setCellInfos(indices_slv, indices_mas, indices_do_slv, indices_do_mas)
        self.assertTrue(obj._flag_CellsInfos, "_flag_CellsInfos should be True")

        # - Check cell information (from setCellInfos)
        self.assertEqual(
            obj._indicesSlaveInterface,
            indices_slv,
            "_indicesSlaveInterface attribute was not initialized correctly.",
        )
        self.assertEqual(
            obj._indicesMasterInterface,
            indices_mas,
            "_indicesMasterInterface attribute was not initialized correctly.",
        )
        self.assertEqual(
            obj._indicesSlaveDomain,
            indices_do_slv,
            "_indicesSlaveDomain attribute was not initialized correctly.",
        )
        self.assertEqual(
            obj._indicesMasterDomain,
            indices_do_mas,
            "_indicesMasterDomain attribute was not initialized correctly.",
        )


class TestPairingAnalysisAster(unittest.TestCase):

    def test_init_TestPairingAnalysisAster(self):
        """Test initialization with a mock (if aster pairing has been computed)"""
        mock_aster_process = unittest.mock.Mock()
        mock_aster_process._hasRun = True
        mock_aster_process._coords = np.array([[0.0, 0.0]])
        mock_aster_process._asterConnectivity = [[0]]
        mock_aster_process._listPairs = np.array([[0, 1]])
        mock_aster_process._intePointsList = np.array([[0.5, 0.5]])
        mock_aster_process._quadPointsList = np.array([[0.25, 0.25]])

        # - Mock the getCells method
        mock_aster_process._asterMesh.getCells.side_effect = [
            ["10"],  # First call _slvSolid
            ["0"],  # Second call _slvtInt
            ["11"],  # Third call for _mastSolid
            ["1"],  # Fourth call for _mastInt
        ]

        obj = PairingAnalysisAster(2, "m_solid", "m_int", "s_solid", "s_int", mock_aster_process)

        self.assertTrue(obj._flag_MeshInfos)
        self.assertTrue(obj._flag_PairingInfos)
        self.assertTrue(obj._flag_CellsInfos)
        np.testing.assert_array_equal(obj._listPairs, mock_aster_process._listPairs)
        self.assertEqual(obj._indicesSlaveDomain, [10])
        self.assertEqual(obj._indicesMasterInterface, [1])

    def test_init_raises_error_if_pairing_not_run(self):
        """Test initialization with a mock (if aster pairing has not been computed)"""
        mock_aster_process = unittest.mock.Mock()
        mock_aster_process._hasRun = False

        with self.assertRaisesRegex(ValueError, "No pairing has been computed before"):
            PairingAnalysisAster(2, "m", "mi", "s", "si", mock_aster_process)


if __name__ == "__main__":
    unittest.main()
