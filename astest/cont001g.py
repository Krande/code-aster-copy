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
from libaster import ContactPairing, CoordinatesSpace


DEBUT(CODE="OUI", ERREUR=_F(ALARME="ALARME"), DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

meshLine = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

meshLine = MODI_MAILLAGE(
    reuse=meshLine, MAILLAGE=meshLine, ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS"))
)

modelLine = AFFE_MODELE(
    MAILLAGE=meshLine, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN")
)

# Slave side - CONT_BAS
contactLine = DEFI_CONT(
    MODELE=modelLine,
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
pair = ContactPairing(contactLine)
pair.compute()

# One zone
test.assertEqual(pair.getNumberOfZones(), 1)
indexZone = 0

# Get number of pairs on first zone
test.assertEqual(pair.getNumberOfPairs(indexZone), 4)

# Get list of pairs on first zone
test.assertSequenceEqual(pair.getListOfPairs(indexZone), [(18, 3), (18, 4), (19, 4), (19, 5)])

# Get intersection points: 4 paires de points d'intersection
test.assertEqual(len(pair.getNumberOfIntersectionPoints(indexZone)), 4)

# For each pair: two points of intersections
iPair = 0
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair], 2)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 1], 2)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 2], 2)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 3], 2)

# Values of intersections
iPair = 0
test.assertSequenceEqual(
    pair.getIntersectionPoints(indexZone, CoordinatesSpace.Slave)[iPair],
    (0.5999999999999999, 0.0, -0.7333333333333334, 0.0),
)
test.assertSequenceEqual(
    pair.getIntersectionPoints(indexZone, CoordinatesSpace.Slave)[iPair + 1],
    (-0.7333333333333334, 0.0, -1.0, 0.0),
)

meshQuad = CREA_MAILLAGE(MAILLAGE=meshLine, LINE_QUAD=_F(TOUT="OUI"))

modelQuad = AFFE_MODELE(
    MAILLAGE=meshQuad, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN")
)

# Slave side - CONT_BAS
contactQuad = DEFI_CONT(
    MODELE=modelQuad,
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
pair = ContactPairing(contactQuad)
pair.compute()

# One zone
test.assertEqual(pair.getNumberOfZones(), 1)
indexZone = 0

# Get number of pairs on first zone
test.assertEqual(pair.getNumberOfPairs(indexZone), 4)

# Get list of pairs on first zone
test.assertSequenceEqual(pair.getListOfPairs(indexZone), [(18, 3), (18, 4), (19, 4), (19, 5)])

# Get intersection points
test.assertEqual(len(pair.getNumberOfIntersectionPoints(indexZone)), 4)

# For each pair: two points of intersections
iPair = 0
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair], 2)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 1], 2)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 2], 2)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 3], 2)

# Values of intersections
iPair = 0
test.assertSequenceEqual(
    pair.getIntersectionPoints(indexZone, CoordinatesSpace.Slave)[iPair],
    (0.5999999999999999, 0.0, -0.7333333333333334, 0.0),
)
test.assertSequenceEqual(
    pair.getIntersectionPoints(indexZone, CoordinatesSpace.Slave)[iPair + 1],
    (-0.7333333333333334, 0.0, -1.0, 0.0),
)


FIN()
