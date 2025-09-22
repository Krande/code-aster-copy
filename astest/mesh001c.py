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

import os.path as osp


from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import SharedTmpdir

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

# check ParallelMesh object API
test = CA.TestCase()
rank = CA.MPI.ASTER_COMM_WORLD.Get_rank()

# from MED format
mesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

pmesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH", INFO=2)

pmesh = DEFI_GROUP(reuse=pmesh, MAILLAGE=pmesh, CREA_GROUP_NO=(_F(TOUT_GROUP_MA="OUI"),))

test.assertTrue(pmesh.isParallel())
test.assertEqual(pmesh.getDimension(), 3)
test.assertFalse(pmesh.isQuadratic())

# test COPIER on ParallelMesh
pmesh2 = COPIER(CONCEPT=pmesh)
test.assertTrue(pmesh2.isParallel())
test.assertEqual(pmesh2.getDimension(), 3)
test.assertFalse(pmesh2.isQuadratic())

# test dimension
TABG = RECU_TABLE(CO=mesh, NOM_TABLE="CARA_GEOM")
test.assertAlmostEqual(TABG["X_MIN", 1], -1.0)
test.assertAlmostEqual(TABG["X_MAX", 1], 6.0)
test.assertAlmostEqual(TABG["Y_MIN", 1], -6.0)
test.assertAlmostEqual(TABG["Y_MAX", 1], 0.0)
test.assertAlmostEqual(TABG["Z_MIN", 1], -1.0)
test.assertAlmostEqual(TABG["Z_MAX", 1], 6.0)
test.assertAlmostEqual(TABG["AR_MIN", 1], 0.3)
test.assertAlmostEqual(TABG["AR_MAX", 1], 1.0)

TABGP = RECU_TABLE(CO=pmesh, NOM_TABLE="CARA_GEOM")
test.assertAlmostEqual(TABG["X_MIN", 1], TABGP["X_MIN", 1])
test.assertAlmostEqual(TABG["X_MAX", 1], TABGP["X_MAX", 1])
test.assertAlmostEqual(TABG["Y_MIN", 1], TABGP["Y_MIN", 1])
test.assertAlmostEqual(TABG["Y_MAX", 1], TABGP["Y_MAX", 1])
test.assertAlmostEqual(TABG["Z_MIN", 1], TABGP["Z_MIN", 1])
test.assertAlmostEqual(TABG["Z_MAX", 1], TABGP["Z_MAX", 1])
test.assertAlmostEqual(TABG["AR_MIN", 1], TABGP["AR_MIN", 1])
test.assertAlmostEqual(TABG["AR_MAX", 1], TABGP["AR_MAX", 1])

global_grp = mesh.getGroupsOfCells()
test.assertSequenceEqual(sorted(pmesh.getGroupsOfNodes()), sorted(global_grp))
test.assertTrue(pmesh.hasGroupOfNodes("TOUT"))

pmesh.printMedFile("mesh_%d.med" % rank)
with SharedTmpdir("mesh001c_") as tmpdir:
    pmesh.printMedFile(osp.join(tmpdir.path, "mesh001c.med"), local=False)

pmesh2 = CA.ParallelMesh()
pmesh2.readMedFile("mesh_%d.med" % rank, partitioned=True)
with SharedTmpdir("mesh001c_") as tmpdir:
    pmesh2.printMedFile(osp.join(tmpdir.path, "mesh001c.med"), local=False)

# test global numbering of nodes
nodes_gnum = pmesh.getNodes(localNumbering=False)
nodes2_gnum = pmesh2.getNodes(localNumbering=False)
test.assertSequenceEqual(nodes_gnum, nodes2_gnum)

test.printSummary()

CA.close()
