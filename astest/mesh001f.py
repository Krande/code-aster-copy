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

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

rank = CA.MPI.ASTER_COMM_WORLD.Get_rank()

test = CA.TestCase()
mesh = CA.ParallelMesh()
mesh.readMedFile("mesh001f.med")
print(mesh.getCoordinates().getValues())

nb_nodes_lin = [6, 6]
nb_cells = [7, 7]
nodes_rank_lin = [[0, 0, 1, 0, 1, 0], [1, 0, 1, 0, 1, 1]]
cells_rank = [[0, 0, 0, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0, 1]]
nodes_nume = [[2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 0, 1]]
nodes_coord = [
    [1.5, -0.5, 0.0, 1.5, 0.5, 0.0, -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, -0.5, 0.5, 0.0, 0.5, 0.5, 0.0],
    [
        -0.5,
        -0.5,
        0.0,
        0.5,
        -0.5,
        0.0,
        -0.5,
        0.5,
        0.0,
        0.5,
        0.5,
        0.0,
        -1.5,
        -0.5,
        0.0,
        -1.5,
        0.5,
        0.0,
    ],
]
test.assertEqual(nb_nodes_lin[rank], mesh.getNumberOfNodes())
test.assertEqual(nb_cells[rank], mesh.getNumberOfCells())
test.assertSequenceEqual(nodes_rank_lin[rank], mesh.getNodesOwner())
test.assertSequenceEqual(cells_rank[rank], mesh.getCellsOwner())
test.assertSequenceEqual(nodes_nume[rank], mesh.getNodes(localNumbering=False))
test.assertEqual(nodes_coord[rank], mesh.getCoordinates().getValues())
test.assertEqual(nodes_coord[rank], mesh.getCoordinates().toFieldOnNodes(mesh).getValues())
coordf = CREA_CHAMP(OPERATION="EXTR", TYPE_CHAM="NOEU_GEOM_R", NOM_CHAM="GEOMETRIE", MAILLAGE=mesh)
test.assertEqual(nodes_coord[rank], coordf.getValues())

mesh_quad = CREA_MAILLAGE(MAILLAGE=mesh, LINE_QUAD=_F(TOUT="OUI"))
nb_nodes_quad = [13, 13]
nodes_rank_quad = [[0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0]]
nodes_nume_quad = [
    [0, 1, 10, 2, 11, 3, 4, 5, 6, 7, 8, 17, 9],
    [10, 2, 11, 3, 12, 13, 14, 7, 15, 16, 4, 17, 9],
]
test.assertEqual(nb_nodes_quad[rank], mesh_quad.getNumberOfNodes())
test.assertEqual(nb_cells[rank], mesh_quad.getNumberOfCells())
test.assertTrue(mesh_quad.isParallel())
test.assertEqual(mesh_quad.getDimension(), 2)
test.assertSequenceEqual(nodes_rank_quad[rank], mesh_quad.getNodesOwner())
test.assertSequenceEqual(cells_rank[rank], mesh_quad.getCellsOwner())
test.assertSequenceEqual(nodes_nume_quad[rank], mesh_quad.getNodes(localNumbering=False))

mesh_biquad = mesh_quad.convertToBiQuadratic()
mesh_lin = mesh_quad.convertToLinear()
test.assertEqual(nb_nodes_lin[rank], mesh_lin.getNumberOfNodes())
test.assertEqual(nb_cells[rank], mesh_lin.getNumberOfCells())
test.assertTrue(mesh_lin.isParallel())
test.assertEqual(mesh_lin.getDimension(), 2)
test.assertSequenceEqual(nodes_rank_lin[rank], mesh_lin.getNodesOwner())
test.assertSequenceEqual(cells_rank[rank], mesh_lin.getCellsOwner())
nodes_nume_lin = [[0, 1, 4, 2, 5, 3], [4, 2, 5, 3, 6, 7]]
test.assertSequenceEqual(nodes_nume_lin[rank], mesh_lin.getNodes(localNumbering=False))

test.printSummary()

CA.close()
