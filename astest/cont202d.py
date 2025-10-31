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
        "nbPairs": 102,
        "listPairs": [
            (186, 422),
            (186, 424),
            (198, 424),
            (198, 422),
            (198, 440),
            (198, 439),
            (198, 438),
            (198, 421),
            (187, 421),
            (202, 421),
            (202, 438),
            (202, 439),
            (202, 468),
            (202, 440),
            (202, 470),
            (202, 422),
            (202, 452),
            (202, 472),
            (202, 424),
            (202, 453),
            (202, 443),
            (189, 421),
            (189, 436),
            (189, 420),
            (203, 472),
            (203, 470),
            (203, 468),
            (203, 464),
            (203, 438),
            (203, 455),
            (203, 421),
            (203, 451),
            (197, 443),
            (197, 424),
            (197, 453),
            (197, 466),
            (197, 452),
            (197, 465),
            (197, 472),
            (197, 470),
            (188, 420),
            (201, 420),
            (201, 436),
            (201, 421),
            (201, 460),
            (201, 438),
            (201, 449),
            (201, 464),
            (201, 468),
            (201, 451),
            (201, 455),
            (195, 451),
            (195, 455),
            (195, 464),
            (195, 458),
            (195, 468),
            (195, 471),
            (195, 456),
            (195, 470),
            (195, 473),
            (195, 469),
            (195, 457),
            (195, 472),
            (195, 461),
            (196, 472),
            (196, 466),
            (196, 473),
            (196, 465),
            (196, 467),
            (196, 463),
            (196, 462),
            (190, 420),
            (200, 451),
            (200, 449),
            (200, 435),
            (200, 460),
            (200, 436),
            (200, 420),
            (193, 461),
            (193, 457),
            (193, 456),
            (193, 441),
            (193, 458),
            (193, 435),
            (193, 434),
            (193, 455),
            (193, 451),
            (194, 461),
            (194, 469),
            (194, 459),
            (194, 471),
            (194, 462),
            (194, 473),
            (194, 463),
            (194, 472),
            (194, 467),
            (191, 420),
            (199, 435),
            (199, 451),
            (199, 441),
            (199, 434),
            (192, 434),
        ],
    },
    2: {
        "nbPairs": 102,
        "listPairs": [
            (186, 422),
            (186, 424),
            (198, 424),
            (198, 422),
            (198, 440),
            (198, 439),
            (198, 438),
            (198, 421),
            (187, 421),
            (202, 421),
            (202, 438),
            (202, 439),
            (202, 468),
            (202, 440),
            (202, 470),
            (202, 422),
            (202, 452),
            (202, 472),
            (202, 424),
            (202, 453),
            (202, 443),
            (189, 421),
            (189, 436),
            (189, 420),
            (203, 472),
            (203, 470),
            (203, 468),
            (203, 464),
            (203, 438),
            (203, 455),
            (203, 421),
            (203, 451),
            (197, 443),
            (197, 424),
            (197, 453),
            (197, 466),
            (197, 452),
            (197, 465),
            (197, 472),
            (197, 470),
            (188, 420),
            (201, 420),
            (201, 436),
            (201, 421),
            (201, 460),
            (201, 438),
            (201, 449),
            (201, 464),
            (201, 468),
            (201, 451),
            (201, 455),
            (195, 451),
            (195, 455),
            (195, 464),
            (195, 458),
            (195, 468),
            (195, 471),
            (195, 456),
            (195, 470),
            (195, 473),
            (195, 469),
            (195, 457),
            (195, 472),
            (195, 461),
            (196, 472),
            (196, 466),
            (196, 473),
            (196, 465),
            (196, 467),
            (196, 463),
            (196, 462),
            (190, 420),
            (200, 451),
            (200, 449),
            (200, 435),
            (200, 460),
            (200, 436),
            (200, 420),
            (193, 461),
            (193, 457),
            (193, 456),
            (193, 441),
            (193, 458),
            (193, 435),
            (193, 434),
            (193, 455),
            (193, 451),
            (194, 461),
            (194, 469),
            (194, 459),
            (194, 471),
            (194, 462),
            (194, 473),
            (194, 463),
            (194, 472),
            (194, 467),
            (191, 420),
            (199, 435),
            (199, 451),
            (199, 441),
            (199, 434),
            (192, 434),
        ],
    },
    3: {
        "nbPairs": 102,
        "listPairs": [
            (186, 422),
            (186, 424),
            (187, 421),
            (188, 420),
            (189, 420),
            (189, 421),
            (189, 436),
            (190, 420),
            (191, 420),
            (192, 434),
            (193, 434),
            (193, 435),
            (193, 441),
            (193, 451),
            (193, 455),
            (193, 456),
            (193, 457),
            (193, 458),
            (193, 461),
            (194, 459),
            (194, 461),
            (194, 462),
            (194, 463),
            (194, 467),
            (194, 469),
            (194, 471),
            (194, 472),
            (194, 473),
            (195, 451),
            (195, 455),
            (195, 456),
            (195, 457),
            (195, 458),
            (195, 461),
            (195, 464),
            (195, 468),
            (195, 469),
            (195, 470),
            (195, 471),
            (195, 472),
            (195, 473),
            (196, 462),
            (196, 463),
            (196, 465),
            (196, 466),
            (196, 467),
            (196, 472),
            (196, 473),
            (197, 424),
            (197, 443),
            (197, 452),
            (197, 453),
            (197, 465),
            (197, 466),
            (197, 470),
            (197, 472),
            (198, 421),
            (198, 422),
            (198, 424),
            (198, 438),
            (198, 439),
            (198, 440),
            (199, 434),
            (199, 435),
            (199, 441),
            (199, 451),
            (200, 420),
            (200, 435),
            (200, 436),
            (200, 449),
            (200, 451),
            (200, 460),
            (201, 420),
            (201, 421),
            (201, 436),
            (201, 438),
            (201, 449),
            (201, 451),
            (201, 455),
            (201, 460),
            (201, 464),
            (201, 468),
            (202, 421),
            (202, 422),
            (202, 424),
            (202, 438),
            (202, 439),
            (202, 440),
            (202, 443),
            (202, 452),
            (202, 453),
            (202, 468),
            (202, 470),
            (202, 472),
            (203, 421),
            (203, 438),
            (203, 451),
            (203, 455),
            (203, 464),
            (203, 468),
            (203, 470),
            (203, 472),
        ],
    },
    4: {
        "nbPairs": 102,
        "listPairs": [
            (420, 188),
            (420, 190),
            (420, 189),
            (420, 191),
            (420, 201),
            (420, 200),
            (436, 200),
            (436, 201),
            (436, 189),
            (421, 189),
            (421, 187),
            (421, 201),
            (421, 198),
            (421, 203),
            (421, 202),
            (460, 201),
            (460, 200),
            (438, 202),
            (438, 198),
            (438, 203),
            (438, 201),
            (449, 200),
            (449, 201),
            (464, 201),
            (464, 203),
            (464, 195),
            (439, 198),
            (439, 202),
            (468, 201),
            (468, 203),
            (468, 202),
            (468, 195),
            (451, 201),
            (451, 200),
            (451, 203),
            (451, 199),
            (451, 195),
            (451, 193),
            (455, 195),
            (455, 193),
            (455, 203),
            (455, 201),
            (440, 202),
            (440, 198),
            (470, 195),
            (470, 203),
            (470, 202),
            (470, 197),
            (435, 193),
            (435, 199),
            (435, 200),
            (458, 193),
            (458, 195),
            (422, 198),
            (422, 186),
            (422, 202),
            (452, 202),
            (452, 197),
            (472, 197),
            (472, 202),
            (472, 196),
            (472, 203),
            (472, 194),
            (472, 195),
            (441, 199),
            (441, 193),
            (471, 195),
            (471, 194),
            (456, 195),
            (456, 193),
            (424, 202),
            (424, 198),
            (424, 197),
            (424, 186),
            (453, 197),
            (453, 202),
            (466, 196),
            (466, 197),
            (473, 195),
            (473, 194),
            (473, 196),
            (434, 193),
            (434, 199),
            (434, 192),
            (469, 194),
            (469, 195),
            (457, 193),
            (457, 195),
            (443, 197),
            (443, 202),
            (465, 197),
            (465, 196),
            (467, 196),
            (467, 194),
            (459, 194),
            (461, 195),
            (461, 193),
            (461, 194),
            (463, 194),
            (463, 196),
            (462, 194),
            (462, 196),
        ],
    },
    5: {
        "nbPairs": 102,
        "listPairs": [
            (420, 188),
            (420, 190),
            (420, 189),
            (420, 191),
            (420, 201),
            (420, 200),
            (436, 200),
            (436, 201),
            (436, 189),
            (421, 189),
            (421, 187),
            (421, 201),
            (421, 198),
            (421, 203),
            (421, 202),
            (460, 201),
            (460, 200),
            (438, 202),
            (438, 198),
            (438, 203),
            (438, 201),
            (449, 200),
            (449, 201),
            (464, 201),
            (464, 203),
            (464, 195),
            (439, 198),
            (439, 202),
            (468, 201),
            (468, 203),
            (468, 202),
            (468, 195),
            (451, 201),
            (451, 200),
            (451, 203),
            (451, 199),
            (451, 195),
            (451, 193),
            (455, 195),
            (455, 193),
            (455, 203),
            (455, 201),
            (440, 202),
            (440, 198),
            (470, 195),
            (470, 203),
            (470, 202),
            (470, 197),
            (435, 193),
            (435, 199),
            (435, 200),
            (458, 193),
            (458, 195),
            (422, 198),
            (422, 186),
            (422, 202),
            (452, 202),
            (452, 197),
            (472, 197),
            (472, 202),
            (472, 196),
            (472, 203),
            (472, 194),
            (472, 195),
            (441, 199),
            (441, 193),
            (471, 195),
            (471, 194),
            (456, 195),
            (456, 193),
            (424, 202),
            (424, 198),
            (424, 197),
            (424, 186),
            (453, 197),
            (453, 202),
            (466, 196),
            (466, 197),
            (473, 195),
            (473, 194),
            (473, 196),
            (434, 193),
            (434, 199),
            (434, 192),
            (469, 194),
            (469, 195),
            (457, 193),
            (457, 195),
            (443, 197),
            (443, 202),
            (465, 197),
            (465, 196),
            (467, 196),
            (467, 194),
            (459, 194),
            (461, 195),
            (461, 193),
            (461, 194),
            (463, 194),
            (463, 196),
            (462, 194),
            (462, 196),
        ],
    },
    6: {
        "nbPairs": 102,
        "listPairs": [
            (420, 188),
            (420, 189),
            (420, 190),
            (420, 191),
            (420, 200),
            (420, 201),
            (421, 187),
            (421, 189),
            (421, 198),
            (421, 201),
            (421, 202),
            (421, 203),
            (422, 186),
            (422, 198),
            (422, 202),
            (424, 186),
            (424, 197),
            (424, 198),
            (424, 202),
            (434, 192),
            (434, 193),
            (434, 199),
            (435, 193),
            (435, 199),
            (435, 200),
            (436, 189),
            (436, 200),
            (436, 201),
            (438, 198),
            (438, 201),
            (438, 202),
            (438, 203),
            (439, 198),
            (439, 202),
            (440, 198),
            (440, 202),
            (441, 193),
            (441, 199),
            (443, 197),
            (443, 202),
            (449, 200),
            (449, 201),
            (451, 193),
            (451, 195),
            (451, 199),
            (451, 200),
            (451, 201),
            (451, 203),
            (452, 197),
            (452, 202),
            (453, 197),
            (453, 202),
            (455, 193),
            (455, 195),
            (455, 201),
            (455, 203),
            (456, 193),
            (456, 195),
            (457, 193),
            (457, 195),
            (458, 193),
            (458, 195),
            (459, 194),
            (460, 200),
            (460, 201),
            (461, 193),
            (461, 194),
            (461, 195),
            (462, 194),
            (462, 196),
            (463, 194),
            (463, 196),
            (464, 195),
            (464, 201),
            (464, 203),
            (465, 196),
            (465, 197),
            (466, 196),
            (466, 197),
            (467, 194),
            (467, 196),
            (468, 195),
            (468, 201),
            (468, 202),
            (468, 203),
            (469, 194),
            (469, 195),
            (470, 195),
            (470, 197),
            (470, 202),
            (470, 203),
            (471, 194),
            (471, 195),
            (472, 194),
            (472, 195),
            (472, 196),
            (472, 197),
            (472, 202),
            (472, 203),
            (473, 194),
            (473, 195),
            (473, 196),
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
