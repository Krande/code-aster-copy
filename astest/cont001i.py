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


DEBUT(CODE="OUI", ERREUR=_F(ALARME="ALARME"), DEBUG=_F(SDVERI="NON"), INFO=1)

test = CA.TestCase()

meshLine = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

meshLine = MODI_MAILLAGE(
    reuse=meshLine, MAILLAGE=meshLine, ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS"))
)

modelLine = AFFE_MODELE(
    MAILLAGE=meshLine, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D")
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
pair.setVerbosity(1)
pair.compute()

# One zone
test.assertEqual(pair.getNumberOfZones(), 1)
indexZone = 0

# Get number of pairs on first zone
test.assertEqual(pair.getNumberOfPairs(indexZone), 9)

# Get list of pairs on first zone
test.assertSequenceEqual(
    pair.getListOfPairs(indexZone),
    [(68, 88), (68, 90), (70, 88), (69, 90), (69, 91), (69, 88), (69, 89), (71, 88), (71, 89)],
)

# Get intersection points: 9 paires de points d'intersection
test.assertEqual(len(pair.getNumberOfIntersectionPoints(indexZone)), 9)

# For each pair: 4 points of intersections
iPair = 0
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair], 4)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 1], 4)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 2], 4)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 3], 4)

# Values of intersections
iPair = 0
test.assertSequenceEqual(
    pair.getIntersectionPoints(indexZone, CoordinatesSpace.Slave)[iPair],
    (1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0),
)
iPair = 1
test.assertSequenceEqual(
    pair.getIntersectionPoints(indexZone, CoordinatesSpace.Slave)[iPair],
    (0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 0.0),
)

meshQuad = CREA_MAILLAGE(MAILLAGE=meshLine, LINE_QUAD=_F(TOUT="OUI"))

modelQuad = AFFE_MODELE(
    MAILLAGE=meshQuad, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D")
)

# Slave side - CONT_BAS
contactQuad = DEFI_CONT(
    MODELE=modelQuad,
    INFO=2,
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

# Pairing checks
pair = ContactPairing(contactQuad)
pair.setVerbosity(1)
pair.compute()

# One zone
test.assertEqual(pair.getNumberOfZones(), 1)
indexZone = 0

# Get number of pairs on first zone
test.assertEqual(pair.getNumberOfPairs(indexZone), 9)

# Get list of pairs on first zone
test.assertSequenceEqual(
    pair.getListOfPairs(indexZone),
    [(88, 68), (88, 70), (88, 69), (88, 71), (90, 69), (90, 68), (89, 71), (89, 69), (91, 69)],
)

# Get intersection points: 4 paires de points d'intersection
test.assertEqual(len(pair.getNumberOfIntersectionPoints(indexZone)), 9)

# For each pair: two points of intersections
iPair = 0
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair], 4)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 1], 4)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 2], 4)
test.assertEqual(pair.getNumberOfIntersectionPoints(indexZone)[iPair + 3], 4)

# Values of intersections
iPair = 0
test.assertSequenceEqual(
    pair.getIntersectionPoints(indexZone, CoordinatesSpace.Slave)[iPair],
    (1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0),
)

FIN()
