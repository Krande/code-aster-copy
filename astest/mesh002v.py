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

import os
from collections import Counter

import code_aster
from code_aster.Commands import *
from code_aster.Utilities import logger, shared_tmpdir
from code_aster.Utilities.MedUtils.MEDPartitioner import MEDPartitioner

MPI = code_aster.MPI

code_aster.init("--test")

test = code_aster.TestCase()

rank = code_aster.MPI.COMM_WORLD.Get_rank()

with shared_tmpdir("mesh002v_") as tmpdir:
    logger.info(f"\n--- Common temporary directory for the testcase: {tmpdir}")

    logger.info("splitting the mesh...")
    ms = MEDPartitioner("sdlx400b.mmed")
    ms.partitionMesh(verbose=True)

    # write the mesh in tmpdir
    logger.info(f"writing mesh...")
    ms.writeMesh(tmpdir)

    # add PO1 (need to load sequential mesh)
    logger.info("adding POI1 elements...")
    ms.addPO1(verbose=True)

    pMesh2 = code_aster.ParallelMesh()
    pMesh2.readMedFile(os.path.join(tmpdir, f"sdlx400b_new_{rank}.med"), True)


# STD Mesh for comparaison
mesh = code_aster.Mesh()
mesh.readMedFile("sdlx400b.mmed")


nbNodes = [682,672]
nbCells = [93,160]

test.assertEqual(pMesh2.getNumberOfNodes(), nbNodes[rank])
test.assertEqual(pMesh2.getNumberOfCells(), nbCells[rank])


# tests
group_no_std = mesh.getGroupsOfNodes(local=False)
group_no_gl  = pMesh2.getGroupsOfNodes(local=False)
test.assertSequenceEqual(sorted(group_no_std), sorted(group_no_gl))

group_ma_std = mesh.getGroupsOfCells(local=False)
group_ma_gl  = pMesh2.getGroupsOfCells(local=False)
test.assertSequenceEqual(sorted(group_ma_std), sorted(group_ma_gl))

nb_nodes_std = mesh.getNumberOfNodes()
nb_nodes_lc = len(pMesh2.getInnerNodes())
nb_nodes_gl = MPI.COMM_WORLD.allreduce(nb_nodes_lc, MPI.SUM)
test.assertEqual(nb_nodes_std, nb_nodes_gl)

nb_cells_std = mesh.getNumberOfCells()
cells_rank = pMesh2.getCellsRank()
nb_cells_lc = Counter(cells_rank)[rank]
nb_cells_gl = MPI.COMM_WORLD.allreduce(nb_cells_lc, MPI.SUM)
test.assertEqual(nb_cells_std, nb_cells_gl)


test.printSummary()


FIN()
