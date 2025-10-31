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
from libaster import PairingMethod

# -------------------------------------
#       - Reference solution
# -------------------------------------
"""
    1: SLAVE = CONT_BAS + PairingMethod.Fast
"""
refValues = {
    1: {
        "nbPairs": 36,
        "listPairs": [
            (141, 260),
            (141, 265),
            (144, 260),
            (144, 255),
            (144, 250),
            (142, 265),
            (142, 266),
            (142, 260),
            (142, 267),
            (142, 261),
            (142, 262),
            (147, 250),
            (145, 250),
            (145, 255),
            (145, 251),
            (145, 260),
            (145, 256),
            (145, 252),
            (145, 261),
            (145, 257),
            (145, 262),
            (143, 262),
            (143, 267),
            (143, 263),
            (143, 268),
            (148, 250),
            (148, 251),
            (148, 252),
            (146, 262),
            (146, 263),
            (146, 257),
            (146, 258),
            (146, 252),
            (146, 253),
            (149, 252),
            (149, 253),
        ],
    },
    2: {
        "nbPairs": 36,
        "listPairs": [
            (141, 260),
            (141, 265),
            (144, 260),
            (144, 255),
            (144, 250),
            (142, 265),
            (142, 266),
            (142, 260),
            (142, 267),
            (142, 261),
            (142, 262),
            (147, 250),
            (145, 250),
            (145, 255),
            (145, 251),
            (145, 260),
            (145, 256),
            (145, 252),
            (145, 261),
            (145, 257),
            (145, 262),
            (143, 262),
            (143, 267),
            (143, 263),
            (143, 268),
            (148, 250),
            (148, 251),
            (148, 252),
            (146, 262),
            (146, 263),
            (146, 257),
            (146, 258),
            (146, 252),
            (146, 253),
            (149, 252),
            (149, 253),
        ],
    },
    3: {
        "nbPairs": 36,
        "listPairs": [
            (141, 260),
            (141, 265),
            (142, 260),
            (142, 261),
            (142, 262),
            (142, 265),
            (142, 266),
            (142, 267),
            (143, 262),
            (143, 263),
            (143, 267),
            (143, 268),
            (144, 250),
            (144, 255),
            (144, 260),
            (145, 250),
            (145, 251),
            (145, 252),
            (145, 255),
            (145, 256),
            (145, 257),
            (145, 260),
            (145, 261),
            (145, 262),
            (146, 252),
            (146, 253),
            (146, 257),
            (146, 258),
            (146, 262),
            (146, 263),
            (147, 250),
            (148, 250),
            (148, 251),
            (148, 252),
            (149, 252),
            (149, 253),
        ],
    },
    4: {
        "nbPairs": 36,
        "listPairs": [
            (250, 144),
            (250, 147),
            (250, 145),
            (250, 148),
            (255, 145),
            (255, 144),
            (251, 148),
            (251, 145),
            (260, 144),
            (260, 145),
            (260, 141),
            (260, 142),
            (256, 145),
            (252, 145),
            (252, 148),
            (252, 146),
            (252, 149),
            (265, 142),
            (265, 141),
            (261, 142),
            (261, 145),
            (257, 145),
            (257, 146),
            (253, 149),
            (253, 146),
            (266, 142),
            (262, 145),
            (262, 146),
            (262, 142),
            (262, 143),
            (258, 146),
            (267, 142),
            (267, 143),
            (263, 143),
            (263, 146),
            (268, 143),
        ],
    },
    5: {
        "nbPairs": 36,
        "listPairs": [
            (250, 144),
            (250, 147),
            (250, 145),
            (250, 148),
            (255, 145),
            (255, 144),
            (251, 148),
            (251, 145),
            (260, 144),
            (260, 145),
            (260, 141),
            (260, 142),
            (256, 145),
            (252, 145),
            (252, 148),
            (252, 146),
            (252, 149),
            (265, 142),
            (265, 141),
            (261, 142),
            (261, 145),
            (257, 145),
            (257, 146),
            (253, 149),
            (253, 146),
            (266, 142),
            (262, 145),
            (262, 146),
            (262, 142),
            (262, 143),
            (258, 146),
            (267, 142),
            (267, 143),
            (263, 143),
            (263, 146),
            (268, 143),
        ],
    },
    6: {
        "nbPairs": 36,
        "listPairs": [
            (250, 144),
            (250, 145),
            (250, 147),
            (250, 148),
            (251, 145),
            (251, 148),
            (252, 145),
            (252, 146),
            (252, 148),
            (252, 149),
            (253, 146),
            (253, 149),
            (255, 144),
            (255, 145),
            (256, 145),
            (257, 145),
            (257, 146),
            (258, 146),
            (260, 141),
            (260, 142),
            (260, 144),
            (260, 145),
            (261, 142),
            (261, 145),
            (262, 142),
            (262, 143),
            (262, 145),
            (262, 146),
            (263, 143),
            (263, 146),
            (265, 141),
            (265, 142),
            (266, 142),
            (267, 142),
            (267, 143),
            (268, 143),
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
