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
from code_aster.CA import MPI

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

pMesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH")

model = AFFE_MODELE(MAILLAGE=pMesh, AFFE=_F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"))

rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()

if nbproc == 2:
    localGrpNodes = [tuple("SOURCE"), tuple()]
    localGrpCells = [tuple("ENTREE"), tuple("SORTIE")]

test.assertEqual(pMesh.getDimension(), 3)
test.assertEqual(pMesh.getNumberOfNodes(), 261)
test.assertEqual(pMesh.getNumberOfCells(), 36)
test.assertTrue(pMesh.isParallel())
test.assertTrue(pMesh.getGroupsOfNodes(False), tuple("SOURCE"))
test.assertTrue(pMesh.getGroupsOfCells(False), sorted(["ENTREE", "SORTIE"]))
test.assertTrue(pMesh.getGroupsOfNodes(), localGrpNodes[rank])
test.assertTrue(pMesh.getGroupsOfCells(), localGrpCells[rank])

test.printSummary()


FIN()
