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

import code_aster

from code_aster.Commands import *

code_aster.init("--test")

rank = code_aster.MPI.ASTER_COMM_WORLD.Get_rank()

test = code_aster.TestCase()
mesh = code_aster.ParallelMesh()
mesh.readMedFile("mesh001f.med")
print(mesh.getCoordinates().getValues())

nb_nodes_lin = [6, 6]
nb_cells = [7, 7]
nodes_rank_lin = [[0, 0, 1, 0, 1, 0], [1, 1, 1, 0, 1, 0]]
cells_rank = [[0, 0, 0, 0, 0, 0, 0], [1, 1, 0, 1, 0, 1, 0]]
nodes_nume = [[2, 3, 4, 5, 6, 7], [0, 1, 4, 5, 6, 7]]
test.assertEqual(nb_nodes_lin[rank], mesh.getNumberOfNodes())
test.assertEqual(nb_cells[rank], mesh.getNumberOfCells())
test.assertSequenceEqual(nodes_rank_lin[rank], mesh.getNodesRank())
test.assertSequenceEqual(cells_rank[rank], mesh.getCellsRank())
test.assertSequenceEqual(nodes_nume[rank], mesh.getNodes(localNumbering=False))

mesh_quad = CREA_MAILLAGE(MAILLAGE=mesh, LINE_QUAD=_F(TOUT="OUI"))
nb_nodes_quad = [13, 13]
nodes_rank_quad = [[0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0], [1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0]]
nodes_nume_quad = [
    [0, 1, 12, 2, 13, 3, 4, 5, 6, 7, 8, 17, 9],
    [10, 11, 12, 2, 13, 3, 14, 15, 4, 16, 7, 17, 9],
]
test.assertEqual(nb_nodes_quad[rank], mesh_quad.getNumberOfNodes())
test.assertEqual(nb_cells[rank], mesh_quad.getNumberOfCells())
test.assertTrue(mesh_quad.isParallel())
test.assertEqual(mesh_quad.getDimension(), 2)
test.assertSequenceEqual(nodes_rank_quad[rank], mesh_quad.getNodesRank())
test.assertSequenceEqual(cells_rank[rank], mesh_quad.getCellsRank())
test.assertSequenceEqual(nodes_nume_quad[rank], mesh_quad.getNodes(localNumbering=False))

mesh_lin = CREA_MAILLAGE(MAILLAGE=mesh_quad, QUAD_LINE=_F(TOUT="OUI"))
test.assertEqual(nb_nodes_lin[rank], mesh_lin.getNumberOfNodes())
test.assertEqual(nb_cells[rank], mesh_lin.getNumberOfCells())
test.assertTrue(mesh_lin.isParallel())
test.assertEqual(mesh_lin.getDimension(), 2)
test.assertSequenceEqual(nodes_rank_lin[rank], mesh_lin.getNodesRank())
test.assertSequenceEqual(cells_rank[rank], mesh_lin.getCellsRank())
nodes_nume_lin = [[0, 1, 6, 2, 7, 3], [4, 5, 6, 2, 7, 3]]
test.assertSequenceEqual(nodes_nume_lin[rank], mesh_lin.getNodes(localNumbering=False))

test.printSummary()

code_aster.close()
