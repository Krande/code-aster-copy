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
from enum import Enum
from code_aster.MacroCommands.defi_cont import DEFI_CONT
from libaster import ContactPairing, ContactComputation, PairingMethod
import numpy

# -------------------------------------
#       - Reference solution
# -------------------------------------
"""
    1: SLAVE = CONT_BAS + PairingMethod.Fast
"""
refValues = {
    1: {
        "nbPairs": 28,
        "listPairs": [
            (12, 77),
            (12, 84),
            (12, 102),
            (47, 102),
            (47, 106),
            (47, 115),
            (47, 107),
            (47, 113),
            (47, 108),
            (47, 111),
            (49, 111),
            (49, 108),
            (49, 113),
            (49, 112),
            (49, 109),
            (49, 114),
            (49, 116),
            (22, 116),
            (22, 114),
            (22, 117),
            (22, 109),
            (22, 121),
            (22, 110),
            (22, 130),
            (19, 135),
            (19, 105),
            (19, 101),
            (19, 89),
        ],
    },
    2: {
        "nbPairs": 28,
        "listPairs": [
            (12, 77),
            (12, 84),
            (12, 102),
            (47, 102),
            (47, 106),
            (47, 115),
            (47, 107),
            (47, 113),
            (47, 108),
            (47, 111),
            (49, 111),
            (49, 108),
            (49, 113),
            (49, 112),
            (49, 109),
            (49, 114),
            (49, 116),
            (22, 116),
            (22, 114),
            (22, 117),
            (22, 109),
            (22, 121),
            (22, 110),
            (22, 130),
            (19, 135),
            (19, 105),
            (19, 101),
            (19, 89),
        ],
    },
    3: {
        "nbPairs": 31,
        "listPairs": [
            (12, 77),
            (12, 84),
            (12, 102),
            (13, 76),
            (13, 83),
            (13, 101),
            (19, 89),
            (19, 101),
            (19, 105),
            (19, 135),
            (22, 109),
            (22, 110),
            (22, 114),
            (22, 116),
            (22, 117),
            (22, 121),
            (22, 130),
            (47, 102),
            (47, 106),
            (47, 107),
            (47, 108),
            (47, 111),
            (47, 113),
            (47, 115),
            (49, 108),
            (49, 109),
            (49, 111),
            (49, 112),
            (49, 113),
            (49, 114),
            (49, 116),
        ],
    },
    4: {
        "nbPairs": 26,
        "listPairs": [
            (76, 13),
            (76, 20),
            (89, 25),
            (89, 29),
            (89, 35),
            (89, 19),
            (110, 22),
            (110, 64),
            (110, 63),
            (110, 65),
            (109, 63),
            (109, 62),
            (109, 64),
            (109, 58),
            (109, 22),
            (109, 49),
            (108, 49),
            (108, 47),
            (108, 58),
            (108, 38),
            (108, 43),
            (108, 57),
            (107, 57),
            (107, 38),
            (107, 34),
            (107, 47),
        ],
    },
    5: {
        "nbPairs": 26,
        "listPairs": [
            (76, 13),
            (76, 20),
            (89, 25),
            (89, 29),
            (89, 35),
            (89, 19),
            (110, 22),
            (110, 64),
            (110, 63),
            (110, 65),
            (109, 63),
            (109, 62),
            (109, 64),
            (109, 58),
            (109, 22),
            (109, 49),
            (108, 49),
            (108, 47),
            (108, 58),
            (108, 38),
            (108, 43),
            (108, 57),
            (107, 57),
            (107, 38),
            (107, 34),
            (107, 47),
        ],
    },
    6: {
        "nbPairs": 28,
        "listPairs": [
            (76, 13),
            (76, 20),
            (77, 12),
            (77, 18),
            (89, 19),
            (89, 25),
            (89, 29),
            (89, 35),
            (107, 34),
            (107, 38),
            (107, 47),
            (107, 57),
            (108, 38),
            (108, 43),
            (108, 47),
            (108, 49),
            (108, 57),
            (108, 58),
            (109, 22),
            (109, 49),
            (109, 58),
            (109, 62),
            (109, 63),
            (109, 64),
            (110, 22),
            (110, 63),
            (110, 64),
            (110, 65),
        ],
    },
}
# -------------------------------------
#       - Read mesh and set model
# -------------------------------------
DEBUT(CODE="OUI", ERREUR=_F(ALARME="ALARME"), INFO=1)

test = CA.TestCase()

Mail = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

Mail = MODI_MAILLAGE(
    reuse=Mail, MAILLAGE=Mail, ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS"))
)

MODI = AFFE_MODELE(MAILLAGE=Mail, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

# -------------------------------------
#       - CONFIG: SLAVE="CONT_BAS" + PairingMethod.Fast
# -------------------------------------
# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.Fast)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
test.assertEqual(refValues[1]["nbPairs"], nbPairs)

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    test.assertSequenceEqual(listPairs, refValues[2]["listPairs"])

# -------------------------------------
#       - CONFIG: SLAVE="CONT_BAS" + PairingMethod.Legacy
# -------------------------------------
# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.Legacy)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
test.assertEqual(refValues[2]["nbPairs"], nbPairs)

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    test.assertSequenceEqual(listPairs, refValues[2]["listPairs"])

# -------------------------------------
#       - CONFIG: SLAVE="CONT_BAS" + PairingMethod.BrutForce
# -------------------------------------
# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.BrutForce)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
test.assertEqual(refValues[3]["nbPairs"], nbPairs)

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    test.assertSequenceEqual(listPairs, refValues[3]["listPairs"])

# -------------------------------------
#       - CONFIG: SLAVE="CONT_HAUT" + PairingMethod.Fast
# -------------------------------------
# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_HAUT", "CONT_BAS")
meshPair.setMethod(PairingMethod.Fast)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
test.assertEqual(refValues[4]["nbPairs"], nbPairs)

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    test.assertSequenceEqual(listPairs, refValues[4]["listPairs"])

# -------------------------------------
#       - CONFIG: SLAVE="CONT_HAUT" + PairingMethod.Legacy
# -------------------------------------
# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_HAUT", "CONT_BAS")
meshPair.setMethod(PairingMethod.Legacy)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
test.assertEqual(refValues[5]["nbPairs"], nbPairs)

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    test.assertSequenceEqual(listPairs, refValues[5]["listPairs"])

# -------------------------------------
#       - CONFIG: SLAVE="CONT_HAUT" + PairingMethod.BrutForce
# -------------------------------------
# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_HAUT", "CONT_BAS")
meshPair.setMethod(PairingMethod.BrutForce)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
test.assertEqual(refValues[6]["nbPairs"], nbPairs)

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    test.assertSequenceEqual(listPairs, refValues[6]["listPairs"])

FIN()
