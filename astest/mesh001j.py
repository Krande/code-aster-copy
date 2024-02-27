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
test.assertEqual(nb_nodes_lin[rank], mesh.getNumberOfNodes())
test.assertEqual(nb_cells[rank], mesh.getNumberOfCells())
test.assertSequenceEqual(nodes_rank_lin[rank], mesh.getNodesOwner())
test.assertSequenceEqual(cells_rank[rank], mesh.getCellsOwner())
test.assertSequenceEqual(nodes_nume[rank], mesh.getNodes(localNumbering=False))

mesh_raf = CREA_MAILLAGE(MAILLAGE=mesh, RAFFINEMENT=_F(TOUT="OUI"))
nb_nodes_raf = [15, 12]
nb_cells_raf = [18, 14]
nodes_rank_raf = [
    [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1],
]
nodes_nume_raf = [
    [0, 1, 12, 2, 13, 3, 4, 5, 6, 7, 8, 19, 9, 10, 11],
    [12, 13, 14, 15, 16, 7, 17, 18, 4, 19, 10, 20],
]
cells_rank_raf = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1],
]
test.assertEqual(nb_nodes_raf[rank], mesh_raf.getNumberOfNodes())
test.assertEqual(nb_cells_raf[rank], mesh_raf.getNumberOfCells())
test.assertTrue(mesh_raf.isParallel())
test.assertEqual(mesh_raf.getDimension(), 2)
test.assertSequenceEqual(nodes_rank_raf[rank], mesh_raf.getNodesOwner())
test.assertSequenceEqual(cells_rank_raf[rank], mesh_raf.getCellsOwner())
test.assertSequenceEqual(nodes_nume_raf[rank], mesh_raf.getNodes(localNumbering=False))

test.printSummary()

CA.close()
