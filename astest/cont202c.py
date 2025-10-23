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
        "nbPairs": 57,
        "listPairs": [
            (186, 319),
            (198, 319),
            (198, 314),
            (198, 309),
            (187, 309),
            (202, 309),
            (202, 314),
            (202, 310),
            (202, 319),
            (202, 315),
            (202, 320),
            (202, 316),
            (189, 309),
            (189, 304),
            (203, 316),
            (203, 315),
            (203, 311),
            (203, 310),
            (203, 306),
            (203, 309),
            (203, 305),
            (197, 316),
            (197, 315),
            (197, 321),
            (197, 320),
            (197, 319),
            (188, 304),
            (201, 304),
            (201, 309),
            (201, 305),
            (201, 310),
            (201, 306),
            (195, 306),
            (195, 311),
            (195, 316),
            (195, 312),
            (195, 317),
            (196, 321),
            (196, 322),
            (196, 316),
            (190, 304),
            (200, 306),
            (200, 305),
            (200, 304),
            (193, 317),
            (193, 312),
            (193, 311),
            (193, 307),
            (193, 306),
            (194, 317),
            (194, 316),
            (194, 322),
            (194, 321),
            (191, 304),
            (199, 306),
            (199, 307),
            (192, 307),
        ],
    },
    2: {
        "nbPairs": 57,
        "listPairs": [
            (186, 319),
            (198, 319),
            (198, 314),
            (198, 309),
            (187, 309),
            (202, 309),
            (202, 314),
            (202, 310),
            (202, 319),
            (202, 315),
            (202, 320),
            (202, 316),
            (189, 309),
            (189, 304),
            (203, 316),
            (203, 315),
            (203, 311),
            (203, 310),
            (203, 306),
            (203, 309),
            (203, 305),
            (197, 316),
            (197, 315),
            (197, 321),
            (197, 320),
            (197, 319),
            (188, 304),
            (201, 304),
            (201, 309),
            (201, 305),
            (201, 310),
            (201, 306),
            (195, 306),
            (195, 311),
            (195, 316),
            (195, 312),
            (195, 317),
            (196, 321),
            (196, 322),
            (196, 316),
            (190, 304),
            (200, 306),
            (200, 305),
            (200, 304),
            (193, 317),
            (193, 312),
            (193, 311),
            (193, 307),
            (193, 306),
            (194, 317),
            (194, 316),
            (194, 322),
            (194, 321),
            (191, 304),
            (199, 306),
            (199, 307),
            (192, 307),
        ],
    },
    3: {
        "nbPairs": 57,
        "listPairs": [
            (186, 319),
            (187, 309),
            (188, 304),
            (189, 304),
            (189, 309),
            (190, 304),
            (191, 304),
            (192, 307),
            (193, 306),
            (193, 307),
            (193, 311),
            (193, 312),
            (193, 317),
            (194, 316),
            (194, 317),
            (194, 321),
            (194, 322),
            (195, 306),
            (195, 311),
            (195, 312),
            (195, 316),
            (195, 317),
            (196, 316),
            (196, 321),
            (196, 322),
            (197, 315),
            (197, 316),
            (197, 319),
            (197, 320),
            (197, 321),
            (198, 309),
            (198, 314),
            (198, 319),
            (199, 306),
            (199, 307),
            (200, 304),
            (200, 305),
            (200, 306),
            (201, 304),
            (201, 305),
            (201, 306),
            (201, 309),
            (201, 310),
            (202, 309),
            (202, 310),
            (202, 314),
            (202, 315),
            (202, 316),
            (202, 319),
            (202, 320),
            (203, 305),
            (203, 306),
            (203, 309),
            (203, 310),
            (203, 311),
            (203, 315),
            (203, 316),
        ],
    },
    4: {
        "nbPairs": 57,
        "listPairs": [
            (304, 188),
            (304, 190),
            (304, 189),
            (304, 191),
            (304, 201),
            (304, 200),
            (309, 201),
            (309, 189),
            (309, 203),
            (309, 187),
            (309, 202),
            (309, 198),
            (305, 200),
            (305, 201),
            (305, 203),
            (314, 198),
            (314, 202),
            (310, 202),
            (310, 203),
            (310, 201),
            (306, 203),
            (306, 201),
            (306, 195),
            (306, 200),
            (306, 193),
            (306, 199),
            (319, 202),
            (319, 198),
            (319, 197),
            (319, 186),
            (315, 202),
            (315, 203),
            (315, 197),
            (311, 203),
            (311, 195),
            (311, 193),
            (307, 199),
            (307, 192),
            (307, 193),
            (320, 197),
            (320, 202),
            (316, 197),
            (316, 202),
            (316, 196),
            (316, 203),
            (316, 194),
            (316, 195),
            (312, 193),
            (312, 195),
            (321, 197),
            (321, 196),
            (321, 194),
            (317, 195),
            (317, 193),
            (317, 194),
            (322, 194),
            (322, 196),
        ],
    },
    5: {
        "nbPairs": 57,
        "listPairs": [
            (304, 188),
            (304, 190),
            (304, 189),
            (304, 191),
            (304, 201),
            (304, 200),
            (309, 201),
            (309, 189),
            (309, 203),
            (309, 187),
            (309, 202),
            (309, 198),
            (305, 200),
            (305, 201),
            (305, 203),
            (314, 198),
            (314, 202),
            (310, 202),
            (310, 203),
            (310, 201),
            (306, 203),
            (306, 201),
            (306, 195),
            (306, 200),
            (306, 193),
            (306, 199),
            (319, 202),
            (319, 198),
            (319, 197),
            (319, 186),
            (315, 202),
            (315, 203),
            (315, 197),
            (311, 203),
            (311, 195),
            (311, 193),
            (307, 199),
            (307, 192),
            (307, 193),
            (320, 197),
            (320, 202),
            (316, 197),
            (316, 202),
            (316, 196),
            (316, 203),
            (316, 194),
            (316, 195),
            (312, 193),
            (312, 195),
            (321, 197),
            (321, 196),
            (321, 194),
            (317, 195),
            (317, 193),
            (317, 194),
            (322, 194),
            (322, 196),
        ],
    },
    6: {
        "nbPairs": 57,
        "listPairs": [
            (304, 188),
            (304, 189),
            (304, 190),
            (304, 191),
            (304, 200),
            (304, 201),
            (305, 200),
            (305, 201),
            (305, 203),
            (306, 193),
            (306, 195),
            (306, 199),
            (306, 200),
            (306, 201),
            (306, 203),
            (307, 192),
            (307, 193),
            (307, 199),
            (309, 187),
            (309, 189),
            (309, 198),
            (309, 201),
            (309, 202),
            (309, 203),
            (310, 201),
            (310, 202),
            (310, 203),
            (311, 193),
            (311, 195),
            (311, 203),
            (312, 193),
            (312, 195),
            (314, 198),
            (314, 202),
            (315, 197),
            (315, 202),
            (315, 203),
            (316, 194),
            (316, 195),
            (316, 196),
            (316, 197),
            (316, 202),
            (316, 203),
            (317, 193),
            (317, 194),
            (317, 195),
            (319, 186),
            (319, 197),
            (319, 198),
            (319, 202),
            (320, 197),
            (320, 202),
            (321, 194),
            (321, 196),
            (321, 197),
            (322, 194),
            (322, 196),
        ],
    },
}
# -------------------------------------
#       - Read mesh and set model
# -------------------------------------
DEBUT(
    CODE="OUI",
    ERREUR=_F(ALARME="ALARME"),
    # DEBUG=_F(SDVERI='OUI',),
    INFO=1,
)

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
    test.assertSequenceEqual(listPairs, refValues[1]["listPairs"])

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
#       - CONFIG: SLAVE="CONT_BAS" + PairingMethod.Legacy
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
#       - CONFIG: SLAVE="CONT_BAS" + PairingMethod.Legacy
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
#       - CONFIG: SLAVE="CONT_BAS" + PairingMethod.Legacy
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
#       - CONFIG: SLAVE="CONT_BAS" + PairingMethod.Legacy
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
