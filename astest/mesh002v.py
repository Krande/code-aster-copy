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

from collections import Counter


from code_aster.Commands import *
from code_aster import CA

MPI = CA.MPI

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

rank = CA.MPI.ASTER_COMM_WORLD.Get_rank()


pMesh2 = CA.ParallelMesh()
pMesh2.readMedFile("sdlx400b.mmed")


# STD Mesh for comparaison
mesh = CA.Mesh()
mesh.readMedFile("sdlx400b.mmed")


nbNodes = [672, 682]
nbCells = [160, 93]

test.assertEqual(pMesh2.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(pMesh2.getNumberOfCells(), nbCells[rank])


# tests
group_no_std = mesh.getGroupsOfNodes(local=False)
group_no_gl = pMesh2.getGroupsOfNodes(local=False)
test.assertSequenceEqual(sorted(group_no_std), sorted(group_no_gl))

group_ma_std = mesh.getGroupsOfCells(local=False)
group_ma_gl = pMesh2.getGroupsOfCells(local=False)
test.assertSequenceEqual(sorted(group_ma_std), sorted(group_ma_gl))

nb_nodes_std = mesh.getNumberOfNodes()
nb_nodes_lc = len(pMesh2.getInnerNodes())
nb_nodes_gl = MPI.ASTER_COMM_WORLD.allreduce(nb_nodes_lc, MPI.SUM)
test.assertEqual(nb_nodes_std, nb_nodes_gl)

nb_cells_std = mesh.getNumberOfCells()
cells_rank = pMesh2.getCellsOwner()
nb_cells_lc = Counter(cells_rank)[rank]
nb_cells_gl = MPI.ASTER_COMM_WORLD.allreduce(nb_cells_lc, MPI.SUM)
test.assertEqual(nb_cells_std, nb_cells_gl)


test.printSummary()


FIN()
