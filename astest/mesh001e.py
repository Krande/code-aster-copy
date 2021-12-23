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
from code_aster.Commands import DEFI_GROUP
from code_aster.Commands import LIRE_MAILLAGE

code_aster.init("--test")

test = code_aster.TestCase()

mesh = LIRE_MAILLAGE(UNITE=20,FORMAT='MED',)

mesh = DEFI_GROUP(reuse=mesh,
                  MAILLAGE=mesh,
                  CREA_GROUP_NO=(
                      _F(GROUP_MA=('TOP', 'BOTTOM')),
                      _F(NOM='TEST', INTERSEC=('TOP', 'BOTTOM'))
                  ))

gnodes = mesh.getGroupsOfNodes()

gnodes_ref = ['N2', 'N3', 'N4', 'N6', 'N5', 'N7', 'N1', 'N8', 'TOP', 'BOTTOM']

test.assertSequenceEqual( sorted(gnodes_ref), sorted(gnodes) )

gcells = mesh.getGroupsOfCells()

gcells_ref = ['S13', 'S21', 'S51', 'S26', 'S24', 'S37', 'S34', 'S84', 'S75',
              'S56', 'S68', 'S78', 'BACK', 'BOTTOM', 'RIGHT', 'FRONT', 'LEFT',
              'TOP', 'VOLUME']

test.assertSequenceEqual(gcells_ref, gcells)

test.printSummary()

code_aster.close()
