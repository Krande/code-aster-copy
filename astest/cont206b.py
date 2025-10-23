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
        "nbPairs": 77,
        "listPairs": [
            (20, 88),
            (20, 112),
            (20, 101),
            (20, 134),
            (20, 125),
            (20, 144),
            (20, 89),
            (20, 143),
            (20, 113),
            (20, 145),
            (20, 102),
            (20, 135),
            (20, 199),
            (20, 126),
            (20, 165),
            (20, 182),
            (20, 151),
            (20, 183),
            (27, 100),
            (27, 124),
            (27, 111),
            (27, 161),
            (27, 149),
            (27, 228),
            (27, 160),
            (27, 168),
            (27, 237),
            (27, 204),
            (27, 172),
            (27, 236),
            (27, 203),
            (27, 188),
            (27, 263),
            (27, 267),
            (27, 313),
            (27, 298),
            (27, 268),
            (27, 318),
            (30, 299),
            (30, 246),
            (30, 350),
            (30, 237),
            (30, 247),
            (30, 228),
            (30, 232),
            (30, 308),
            (30, 227),
            (30, 161),
            (30, 234),
            (30, 311),
            (30, 219),
            (30, 124),
            (30, 233),
            (30, 310),
            (30, 213),
            (30, 223),
            (57, 223),
            (57, 215),
            (57, 233),
            (57, 202),
            (57, 216),
            (57, 235),
            (57, 99),
            (57, 200),
            (57, 231),
            (57, 309),
            (57, 193),
            (57, 229),
            (57, 242),
            (57, 317),
            (57, 243),
            (57, 315),
            (57, 304),
            (57, 340),
            (57, 341),
            (57, 433),
            (57, 391),
        ],
    },
    2: {
        "nbPairs": 77,
        "listPairs": [
            (20, 88),
            (20, 112),
            (20, 101),
            (20, 134),
            (20, 125),
            (20, 144),
            (20, 89),
            (20, 143),
            (20, 113),
            (20, 145),
            (20, 102),
            (20, 135),
            (20, 199),
            (20, 126),
            (20, 165),
            (20, 182),
            (20, 151),
            (20, 183),
            (27, 100),
            (27, 124),
            (27, 111),
            (27, 161),
            (27, 149),
            (27, 228),
            (27, 160),
            (27, 168),
            (27, 237),
            (27, 204),
            (27, 172),
            (27, 236),
            (27, 203),
            (27, 188),
            (27, 263),
            (27, 267),
            (27, 313),
            (27, 298),
            (27, 268),
            (27, 318),
            (30, 299),
            (30, 246),
            (30, 350),
            (30, 237),
            (30, 247),
            (30, 228),
            (30, 232),
            (30, 308),
            (30, 227),
            (30, 161),
            (30, 234),
            (30, 311),
            (30, 219),
            (30, 124),
            (30, 233),
            (30, 310),
            (30, 213),
            (30, 223),
            (57, 223),
            (57, 215),
            (57, 233),
            (57, 202),
            (57, 216),
            (57, 235),
            (57, 99),
            (57, 200),
            (57, 231),
            (57, 309),
            (57, 193),
            (57, 229),
            (57, 242),
            (57, 317),
            (57, 243),
            (57, 315),
            (57, 304),
            (57, 340),
            (57, 341),
            (57, 433),
            (57, 391),
        ],
    },
    3: {
        "nbPairs": 118,
        "listPairs": [
            (20, 88),
            (20, 89),
            (20, 101),
            (20, 102),
            (20, 112),
            (20, 113),
            (20, 125),
            (20, 126),
            (20, 134),
            (20, 135),
            (20, 143),
            (20, 144),
            (20, 145),
            (20, 151),
            (20, 165),
            (20, 182),
            (20, 183),
            (20, 199),
            (21, 86),
            (21, 87),
            (21, 98),
            (21, 100),
            (21, 110),
            (21, 111),
            (21, 122),
            (21, 123),
            (21, 132),
            (21, 133),
            (21, 148),
            (21, 149),
            (21, 155),
            (21, 168),
            (21, 169),
            (21, 176),
            (21, 186),
            (21, 187),
            (27, 100),
            (27, 111),
            (27, 124),
            (27, 149),
            (27, 160),
            (27, 161),
            (27, 168),
            (27, 172),
            (27, 188),
            (27, 203),
            (27, 204),
            (27, 228),
            (27, 236),
            (27, 237),
            (27, 263),
            (27, 267),
            (27, 268),
            (27, 298),
            (27, 313),
            (27, 318),
            (30, 124),
            (30, 161),
            (30, 213),
            (30, 219),
            (30, 223),
            (30, 227),
            (30, 228),
            (30, 232),
            (30, 233),
            (30, 234),
            (30, 237),
            (30, 246),
            (30, 247),
            (30, 299),
            (30, 308),
            (30, 310),
            (30, 311),
            (30, 350),
            (55, 88),
            (55, 99),
            (55, 112),
            (55, 134),
            (55, 144),
            (55, 177),
            (55, 193),
            (55, 200),
            (55, 209),
            (55, 229),
            (55, 253),
            (55, 256),
            (55, 258),
            (55, 261),
            (55, 270),
            (55, 276),
            (55, 277),
            (55, 293),
            (55, 300),
            (55, 305),
            (55, 337),
            (55, 352),
            (55, 353),
            (57, 99),
            (57, 193),
            (57, 200),
            (57, 202),
            (57, 215),
            (57, 216),
            (57, 223),
            (57, 229),
            (57, 231),
            (57, 233),
            (57, 235),
            (57, 242),
            (57, 243),
            (57, 304),
            (57, 309),
            (57, 315),
            (57, 317),
            (57, 340),
            (57, 341),
            (57, 391),
            (57, 433),
        ],
    },
    4: {
        "nbPairs": 66,
        "listPairs": [
            (87, 21),
            (87, 28),
            (111, 21),
            (111, 27),
            (100, 27),
            (100, 21),
            (149, 27),
            (149, 21),
            (149, 37),
            (149, 28),
            (149, 33),
            (88, 20),
            (88, 55),
            (124, 27),
            (124, 30),
            (160, 37),
            (160, 27),
            (168, 33),
            (112, 55),
            (112, 20),
            (99, 55),
            (99, 57),
            (213, 30),
            (161, 30),
            (161, 27),
            (161, 37),
            (134, 20),
            (134, 55),
            (134, 26),
            (134, 46),
            (134, 42),
            (193, 57),
            (193, 55),
            (193, 46),
            (202, 57),
            (219, 30),
            (219, 72),
            (228, 37),
            (228, 43),
            (228, 27),
            (228, 73),
            (228, 30),
            (228, 72),
            (177, 42),
            (177, 46),
            (177, 55),
            (200, 46),
            (200, 55),
            (200, 65),
            (200, 57),
            (200, 51),
            (200, 66),
            (215, 57),
            (215, 66),
            (223, 72),
            (223, 30),
            (223, 57),
            (223, 66),
            (227, 72),
            (227, 30),
            (227, 73),
            (216, 66),
            (216, 57),
            (216, 51),
            (233, 70),
            (233, 71),
        ],
    },
    5: {
        "nbPairs": 66,
        "listPairs": [
            (87, 21),
            (87, 28),
            (111, 21),
            (111, 27),
            (100, 27),
            (100, 21),
            (149, 27),
            (149, 21),
            (149, 37),
            (149, 28),
            (149, 33),
            (88, 20),
            (88, 55),
            (124, 27),
            (124, 30),
            (160, 37),
            (160, 27),
            (168, 33),
            (112, 55),
            (112, 20),
            (99, 55),
            (99, 57),
            (213, 30),
            (161, 30),
            (161, 27),
            (161, 37),
            (134, 20),
            (134, 55),
            (134, 26),
            (134, 46),
            (134, 42),
            (193, 57),
            (193, 55),
            (193, 46),
            (202, 57),
            (219, 30),
            (219, 72),
            (228, 37),
            (228, 43),
            (228, 27),
            (228, 73),
            (228, 30),
            (228, 72),
            (177, 42),
            (177, 46),
            (177, 55),
            (200, 46),
            (200, 55),
            (200, 65),
            (200, 57),
            (200, 51),
            (200, 66),
            (215, 57),
            (215, 66),
            (223, 72),
            (223, 30),
            (223, 57),
            (223, 66),
            (227, 72),
            (227, 30),
            (227, 73),
            (216, 66),
            (216, 57),
            (216, 51),
            (233, 70),
            (233, 71),
        ],
    },
    6: {
        "nbPairs": 68,
        "listPairs": [
            (87, 21),
            (87, 28),
            (88, 20),
            (88, 55),
            (99, 55),
            (99, 57),
            (100, 21),
            (100, 27),
            (101, 20),
            (101, 26),
            (111, 21),
            (111, 27),
            (112, 20),
            (112, 55),
            (124, 27),
            (124, 30),
            (134, 20),
            (134, 26),
            (134, 42),
            (134, 46),
            (134, 55),
            (149, 21),
            (149, 27),
            (149, 28),
            (149, 33),
            (149, 37),
            (160, 27),
            (160, 37),
            (161, 27),
            (161, 30),
            (161, 37),
            (168, 33),
            (177, 42),
            (177, 46),
            (177, 55),
            (193, 46),
            (193, 55),
            (193, 57),
            (200, 46),
            (200, 51),
            (200, 55),
            (200, 57),
            (200, 65),
            (200, 66),
            (202, 57),
            (213, 30),
            (215, 57),
            (215, 66),
            (216, 51),
            (216, 57),
            (216, 66),
            (219, 30),
            (219, 72),
            (223, 30),
            (223, 57),
            (223, 66),
            (223, 72),
            (227, 30),
            (227, 72),
            (227, 73),
            (228, 27),
            (228, 30),
            (228, 37),
            (228, 43),
            (228, 72),
            (228, 73),
            (233, 70),
            (233, 71),
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
