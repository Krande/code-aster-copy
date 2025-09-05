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


from code_aster.Commands import *
from code_aster import CA

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

mesh_3d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=20)

mesh_3d_raf = CREA_MAILLAGE(MAILLAGE=mesh_3d, RAFFINEMENT=_F(TOUT="OUI", NIVEAU=1), INFO=1)

test.assertEqual(mesh_3d_raf.getDimension(), 3)
test.assertTrue(mesh_3d_raf.isQuadratic())
test.assertEqual(mesh_3d_raf.getNumberOfNodes(), 149)
test.assertEqual(mesh_3d_raf.getNumberOfCells(), 88)


mesh_2d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=21)

mesh_2d_raf = CREA_MAILLAGE(MAILLAGE=mesh_2d, RAFFINEMENT=_F(TOUT="OUI", NIVEAU=1), INFO=1)

test.assertEqual(mesh_2d_raf.getDimension(), 2)
test.assertTrue(mesh_2d_raf.isQuadratic())
test.assertEqual(mesh_2d_raf.getNumberOfNodes(), 40)
test.assertEqual(mesh_2d_raf.getNumberOfCells(), 20)


test.printSummary()

FIN()
