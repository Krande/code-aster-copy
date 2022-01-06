# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
from code_aster.Commands import LIRE_MAILLAGE, CREA_MAILLAGE, FIN

code_aster.init("--test")

test = code_aster.TestCase()

mesh_init = LIRE_MAILLAGE()

mesh = CREA_MAILLAGE(MAILLAGE=mesh_init,
                     RAFFINEMENT=_F(TOUT='OUI', NIVEAU=2,),
                     INFO=1,)


test.assertEqual(mesh.getDimension(), 3)
test.assertTrue(mesh.isQuadratic())
test.assertEqual(mesh.getNumberOfNodes(), 3218)
test.assertEqual(mesh.getNumberOfCells(), 2144)


test.assertSequenceEqual(sorted(mesh.getGroupsOfCells()), [
                         'HAUT', 'LEFT', 'PART', 'PART_extruded', 'PART_top', 'RIGHT'])
test.assertSequenceEqual(mesh.getCells('LEFT'), [1, 2, 3, 4])
test.assertEqual(len(mesh.getCells('HAUT')), 64)

IMPR_RESU(UNITE=50, FORMAT="MED", RESU=_F(MAILLAGE=mesh))
mesh_read = LIRE_MAILLAGE(UNITE=50)

test.printSummary()

FIN()
