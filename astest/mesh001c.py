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

import code_aster

from code_aster.Commands import LIRE_MAILLAGE, DEFI_GROUP

code_aster.init("--test")

# check ParallelMesh object API
test = code_aster.TestCase()
rank = code_aster.MPI.COMM_WORLD.Get_rank()

# from MED format
mesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

pmesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH", INFO = 2)

pmesh=DEFI_GROUP( reuse=pmesh, MAILLAGE=pmesh,
                  CREA_GROUP_NO=(_F(  TOUT_GROUP_MA='OUI'),),)

test.assertTrue(pmesh.isParallel())
test.assertEqual(pmesh.getDimension(), 3)

global_grp = mesh.getGroupsOfCells()
test.assertSequenceEqual(sorted(pmesh.getGroupsOfNodes()), sorted(global_grp))
test.assertTrue( pmesh.hasGroupOfNodes("TOUT"))

pmesh.printMedFile("mesh_%d.med"%rank)

pmesh2 = code_aster.ParallelMesh()
pmesh2.readMedFile("mesh_%d.med"%rank, True)


#test global numbering of nodes
nodes_gnum = pmesh.getNodes(localNumbering=False)
nodes2_gnum = pmesh2.getNodes(localNumbering=False)
test.assertSequenceEqual(nodes_gnum, nodes2_gnum)

test.printSummary()

code_aster.close()
