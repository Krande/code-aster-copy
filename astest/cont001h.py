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

#

from code_aster.Commands import *
from code_aster import CA
from code_aster.MacroCommands.defi_cont import DEFI_CONT
from libaster import ContactPairing


DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

Mail = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

Mail = MODI_MAILLAGE(
    reuse=Mail, MAILLAGE=Mail, ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS"))
)


MODI = AFFE_MODELE(MAILLAGE=Mail, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"))

# Slave side - CONT_BAS
DEFICO_BAS = DEFI_CONT(
    MODELE=MODI,
    INFO=2,
    ZONE=(
        _F(
            APPARIEMENT="MORTAR",
            GROUP_MA_MAIT="CONT_HAUT",
            GROUP_MA_ESCL="CONT_BAS",
            ALGO_CONT="LAGRANGIEN",
            CONTACT_INIT="OUI",
            COEF_MULT_APPA=-1,
        ),
    ),
)

# Definition checks
zone = DEFICO_BAS.getContactZone(0)
test.assertSequenceEqual(zone.getSlaveNodes(), [7, 10, 57, 58, 59, 174, 175, 176, 177])
test.assertSequenceEqual(zone.getSlaveCells(), [60, 61, 62, 63])

# Pairing
pair = ContactPairing(DEFICO_BAS)
pair.compute()

# Check FED creation
fed = pair.getFiniteElementDescriptor()
nema = fed.getVirtualCellsDescriptor()
nemab = []
for k in range(len(nema)):
    nemab.append(nema[k][:-1])
grel = fed.getListOfGroupsOfElements()
test.assertEqual(len(grel), 2)
test.assertEqual(len(grel[0]), 3)
test.assertSequenceEqual(
    nemab,
    [
        [58, 59, 176, 14, 13, 116],
        [58, 59, 176, 13, 1, 115],
        [59, 60, 177, 14, 13, 116],
        [59, 60, 177, 15, 14, 117],
        [60, 8, 178, 15, 14, 117],
        [60, 8, 178, 2, 15, 118],
        [11],
        [175],
    ],
)

# Slave side - CONT_HAUT
DEFICO_HAUT = DEFI_CONT(
    MODELE=MODI,
    INFO=2,
    ZONE=(
        _F(
            APPARIEMENT="MORTAR",
            GROUP_MA_MAIT="CONT_BAS",
            GROUP_MA_ESCL="CONT_HAUT",
            ALGO_CONT="LAGRANGIEN",
            CONTACT_INIT="OUI",
            COEF_MULT_APPA=10e6,
        ),
    ),
)

# Definition checks
zone = DEFICO_HAUT.getContactZone(0)
test.assertSequenceEqual(zone.getSlaveNodes(), [0, 1, 12, 13, 14, 114, 115, 116, 117])
test.assertSequenceEqual(zone.getSlaveCells(), [0, 1, 2, 3])


# Pairing
pair = ContactPairing(DEFICO_HAUT)
pair.compute()

# Check FED creation
fed = pair.getFiniteElementDescriptor()
nema = fed.getVirtualCellsDescriptor()
nemab = []
for k in range(len(nema)):
    nemab.append(nema[k][:-1])
grel = fed.getListOfGroupsOfElements()
test.assertEqual(len(grel), 2)
test.assertEqual(len(grel[0]), 3)
test.assertSequenceEqual(
    nemab,
    [
        [14, 13, 116, 58, 59, 176],
        [14, 13, 116, 11, 58, 175],
        [15, 14, 117, 59, 60, 177],
        [15, 14, 117, 60, 8, 178],
        [2, 15, 118, 60, 8, 178],
        [1],
        [115],
    ],
)


FIN()
