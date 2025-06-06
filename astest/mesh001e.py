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

CA.init("--test")

test = CA.TestCase()

mesh = CA.Mesh.buildCube()

mesh = DEFI_GROUP(
    reuse=mesh,
    MAILLAGE=mesh,
    CREA_GROUP_NO=(_F(GROUP_MA=("TOP", "BOTTOM")), _F(NOM="TEST", INTERSEC=("TOP", "BOTTOM"))),
)

gnodes = sorted(mesh.getGroupsOfNodes())

gnodes_ref = sorted(["N2", "N3", "N4", "N6", "N5", "N7", "N1", "N8", "TOP", "BOTTOM"])

test.assertSequenceEqual(gnodes_ref, gnodes)

gcells = sorted(mesh.getGroupsOfCells())

gcells_ref = sorted(
    [
        "S13",
        "S21",
        "S51",
        "S26",
        "S24",
        "S37",
        "S34",
        "S84",
        "S75",
        "S56",
        "S68",
        "S78",
        "BACK",
        "BOTTOM",
        "RIGHT",
        "FRONT",
        "LEFT",
        "TOP",
        "VOLUME",
    ]
)

test.assertSequenceEqual(gcells_ref, gcells)

test.printSummary()

CA.close()
