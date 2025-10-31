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
## - Pairing test in 2D for linear meshes: TRIA3 vs TRIA3

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
        "nbPairs": 10,
        "listPairs": [
            (2, 34),
            (2, 35),
            (2, 36),
            (2, 37),
            (3, 34),
            (3, 33),
            (3, 32),
            (4, 31),
            (4, 30),
            (4, 29),
        ],
    },
    2: {
        "nbPairs": 10,
        "listPairs": [
            (2, 34),
            (2, 35),
            (2, 36),
            (2, 37),
            (3, 34),
            (3, 33),
            (3, 32),
            (4, 31),
            (4, 30),
            (4, 29),
        ],
    },
    3: {
        "nbPairs": 10,
        "listPairs": [
            (2, 34),
            (2, 35),
            (2, 36),
            (2, 37),
            (3, 32),
            (3, 33),
            (3, 34),
            (4, 29),
            (4, 30),
            (4, 31),
        ],
    },
    4: {"nbPairs": 6, "listPairs": [(31, 4), (31, 5), (32, 4), (32, 3), (33, 3), (34, 2)]},
    5: {"nbPairs": 6, "listPairs": [(31, 4), (31, 5), (32, 4), (32, 3), (33, 3), (34, 2)]},
    6: {"nbPairs": 6, "listPairs": [(31, 4), (31, 5), (32, 3), (32, 4), (33, 3), (34, 2)]},
}
# -------------------------------------
#       - Read mesh and set model
# -------------------------------------
DEBUT(CODE="OUI", ERREUR=_F(ALARME="ALARME"), IGNORE_ALARM=("MODELE1_63", "MODELE1_64"), INFO=1)

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
