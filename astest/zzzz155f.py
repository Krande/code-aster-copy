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


from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA

CA.init("--test")

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
size = MPI.ASTER_COMM_WORLD.Get_size()

pMesh = CA.ParallelMesh()
pMesh.readMedFile("zzzz155f.med")

pMesh = DEFI_GROUP(
    CREA_GROUP_MA=(
        _F(NOM="GRPM1", TYPE_MAILLE="TOUT", UNION=("M_A", "M_B")),
        _F(NOM="GRPM2", TYPE_MAILLE="TOUT", INTERSEC=("M_A U M_B", "M_A")),
        _F(NOM="GRPM3", TYPE_MAILLE="TOUT", GROUP_MA=("M_B",)),
        _F(NOM="GRPM4", TYPE_MAILLE="TOUT", DIFFE=("M_A U M_B", "M_B")),
    ),
    CREA_GROUP_NO=(
        _F(NOM="GRPN1", UNION=("A", "B")),
        _F(NOM="GRPN2", INTERSEC=("A U B", "A")),
        _F(NOM="GRPN3", GROUP_NO=("B",)),
        _F(NOM="GRPN4", DIFFE=("A U B", "B")),
    ),
    MAILLAGE=pMesh,
)

if rank == 0:
    test.assertFalse(pMesh.hasGroupOfNodes("GRPN1", True))
    test.assertFalse(pMesh.hasGroupOfNodes("GRPN2", True))
    test.assertFalse(pMesh.hasGroupOfNodes("GRPN3", True))
    test.assertFalse(pMesh.hasGroupOfNodes("GRPN4", True))
    test.assertFalse(pMesh.hasGroupOfCells("GRPM1", True))
    test.assertFalse(pMesh.hasGroupOfCells("GRPM2", True))
    test.assertFalse(pMesh.hasGroupOfCells("GRPM3", True))
    test.assertFalse(pMesh.hasGroupOfCells("GRPM4", True))
elif rank == 1:
    test.assertTrue(pMesh.hasGroupOfNodes("GRPN1", True))
    test.assertTrue([874, 914] == pMesh.getNodes("GRPN1"))
    test.assertTrue(pMesh.hasGroupOfNodes("GRPN2", True))
    test.assertTrue([874] == pMesh.getNodes("GRPN2"))
    test.assertTrue(pMesh.hasGroupOfNodes("GRPN3", True))
    test.assertTrue([914] == pMesh.getNodes("GRPN3"))
    test.assertTrue(pMesh.hasGroupOfNodes("GRPN4", True))
    test.assertTrue([874] == pMesh.getNodes("GRPN4"))
    test.assertTrue(pMesh.hasGroupOfCells("GRPM1", True))
    test.assertTrue([924, 939], pMesh.getCells("GRPM1"))
    test.assertTrue(pMesh.hasGroupOfCells("GRPM2", True))
    test.assertTrue([924] == pMesh.getCells("GRPM2"))
    test.assertTrue(pMesh.hasGroupOfCells("GRPM3", True))
    test.assertTrue([939] == pMesh.getCells("GRPM3"))
    test.assertTrue(pMesh.hasGroupOfCells("GRPM4", True))
    test.assertTrue([924] == pMesh.getCells("GRPM4"))

test.printSummary()

CA.close()
