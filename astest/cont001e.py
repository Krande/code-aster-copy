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

#
import code_aster
from code_aster.Commands import *
from code_aster.MacroCommands.defi_cont import DEFI_CONT
from libaster import ContactPairing, ContactComputation
import numpy


DEBUT(
    CODE=_F(NIV_PUB_WEB="INTERNET"),
    # DEBUG=_F(SDVERI='OUI',),
    INFO=1,
)

test = code_aster.TestCase()

Mail = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

Mail = MODI_MAILLAGE(
    reuse=Mail,
    MAILLAGE=Mail,
    ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS",)),
)


MODI = AFFE_MODELE(MAILLAGE=Mail,
                   AFFE=_F(TOUT='OUI',
                           PHENOMENE='MECANIQUE',
                           MODELISATION='D_PLAN',),)

# Slave side - CONT_BAS
DEFICO_BAS = DEFI_CONT(
    MODELE=MODI,
    INFO=2,
    LISSAGE="OUI",
    ZONE=(
        _F(
            APPARIEMENT="MORTAR",
            GROUP_MA_MAIT="CONT_HAUT",
            GROUP_MA_ESCL="CONT_BAS",
            ALGO_CONT="LAGRANGIEN",
            CONTACT_INIT="OUI",
        ),
    ),
)

# Dfeinition checks
zone = DEFICO_BAS.getContactZone(0)
test.assertSequenceEqual(zone.getSlaveNodes(), [1, 3, 10, 11])
test.assertSequenceEqual(zone.getSlaveCells(), [9, 10, 11])


# Pairing checks

pair = ContactPairing(DEFICO_BAS)
pair.compute()


# check FED creation
CD = ContactComputation(DEFICO_BAS)
CD.buildContactResFED(pair)
fed = CD.getFiniteElementDescriptor()
nema = fed.getNema()
grel = fed.getListOfGroupOfElements()
test.assertEqual(len(grel), 1)
test.assertEqual(len(grel[0]), 7)
test.assertEqual(len(grel[0]), len(nema) + 1)
test.assertSequenceEqual(nema, [[11, 2, 17, 24, 97],  [11, 2, 24, 25, 97],
                                [12, 11, 24, 25, 97], [12, 11, 25, 26, 97],
                                [12, 11, 26, 19, 97], [4, 12, 26, 19, 97]])

gap, i_gap = CD.geometricGap(pair.getCoordinates())
test.assertEqual(gap.size(), 4)
test.assertEqual(i_gap.size(), 4)
val = gap.getValues()
test.assertTrue(numpy.isnan(val[1]))
val[1] = None
test.assertSequenceEqual(
    val, [0.0, None, 33.33333333333333, 66.66666666666667])
test.assertSequenceEqual(i_gap.getValues(), [1.0, 0.0, 1.0, 1.0])
print(type(gap[1]))

print(gap.getValues())
print(i_gap.getValues())

# gap.printMedFile("/home/C00976/tmp/gap.med")


# Slave side - CONT_HAUT
DEFICO_HAUT = DEFI_CONT(
    MODELE=MODI,
    INFO=2,
    LISSAGE="OUI",
    ZONE=(
        _F(
            APPARIEMENT="MORTAR",
            GROUP_MA_MAIT="CONT_BAS",
            GROUP_MA_ESCL="CONT_HAUT",
            ALGO_CONT="LAGRANGIEN",
            CONTACT_INIT="OUI",
        ),
    ),
)

# Dfeinition checks
zone = DEFICO_HAUT.getContactZone(0)
test.assertSequenceEqual(zone.getSlaveNodes(), [16, 18, 23, 24, 25])
test.assertSequenceEqual(zone.getSlaveCells(), [16, 17, 18, 19])


# Pairing checks

pair = ContactPairing(DEFICO_HAUT)
pair.compute()


# check FED creation
CD = ContactComputation(DEFICO_HAUT)
CD.buildContactResFED(pair)
fed = CD.getFiniteElementDescriptor()
nema = fed.getNema()
grel = fed.getListOfGroupOfElements()
test.assertEqual(len(grel), 2)
test.assertEqual(len(grel[0]), 6)
test.assertEqual(len(grel[1]), 2)
test.assertEqual(len(grel[0]) + len(grel[1]), len(nema) + 2)
test.assertSequenceEqual(nema, [[17, 24, 11, 2, 97],
                                [17, 24, 12, 11, 97],
                                [24, 25, 12, 11, 97],
                                [24, 25, 4, 12, 97],
                                [25, 26, 4, 12, 97],
                                [19, 72]])

gap, i_gap = CD.geometricGap(pair.getCoordinates())
val = gap.getValues()
test.assertTrue(numpy.isnan(val[1]))
val[1] = None
val[4] = None
test.assertEqual(gap.size(), 5)
test.assertEqual(i_gap.size(), 5)
test
test.assertSequenceEqual(
    val, [1.0048591735576161e-14, None, 25.000000000000007, 50.000000000000014, None])
test.assertSequenceEqual(i_gap.getValues(), [1.0, 0.0, 1.0, 1.0, 0.0])

print(gap.getValues())
print(i_gap.getValues())
IMPR_RESU(FORMAT="MED",RESU=(_F(CHAM_GD=gap,)))

data = CD.contactData(pair)
test.assertEqual(data.size() , 60 * len(nema))

FIN()
