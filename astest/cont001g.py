# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
    CODE=_F(NIV_PUB_WEB="INTERNET"), ERREUR=_F(ALARME="EXCEPTION"), DEBUG=_F(SDVERI="OUI"), INFO=1
)

test = code_aster.TestCase()

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
        ),
    ),
)

# Pairing checks

pair = ContactPairing(DEFICO_BAS)
pair.compute()

test.assertSequenceEqual(pair.getListOfPairsOfZone(0), [(18, 3), (18, 4), (19, 4), (19, 5)])
ref = [
    [
        0.5999999999999999,
        -0.7333333333333334,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        -0.7333333333333334,
        -1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        1.0,
        -0.06666666666666687,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        -0.06666666666666687,
        -1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
]

test.assertEqual(len(pair.getSlaveIntersectionPoints(0)), 4)
print(pair.getSlaveIntersectionPoints(0))
test.assertSequenceEqual(pair.getSlaveIntersectionPoints(0), ref)


MailQ = CREA_MAILLAGE(MAILLAGE=Mail, LINE_QUAD=_F(TOUT="OUI"))


MODIQ = AFFE_MODELE(
    MAILLAGE=MailQ, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN")
)

# Slave side - CONT_BAS
DEFICOQ_BAS = DEFI_CONT(
    MODELE=MODIQ,
    INFO=2,
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

# Pairing checks

pair = ContactPairing(DEFICOQ_BAS)
pair.compute()


test.assertSequenceEqual(pair.getListOfPairsOfZone(0), [(18, 3), (18, 4), (19, 4), (19, 5)])
ref = [
    [
        0.5999999999999999,
        -0.7333333333333334,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        -0.7333333333333334,
        -1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        1.0,
        -0.06666666666666687,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    [
        -0.06666666666666687,
        -1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
]

test.assertEqual(len(pair.getSlaveIntersectionPoints(0)), 4)
test.assertEqual(len(pair.getSlaveIntersectionPoints(0)[0]), 16)
print(pair.getSlaveIntersectionPoints(0))
test.assertSequenceEqual(pair.getSlaveIntersectionPoints(0), ref)

FIN()
