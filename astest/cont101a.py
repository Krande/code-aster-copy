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
from libaster import PairingMethod


DEBUT(
    CODE="OUI",
    ERREUR=_F(ALARME="ALARME"),
    # DEBUG=_F(SDVERI='OUI',),
    INFO=1,
)

test = CA.TestCase()

Mail = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", INFO=2)

Mail = MODI_MAILLAGE(
    reuse=Mail, MAILLAGE=Mail, ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS"))
)

MODI = AFFE_MODELE(MAILLAGE=Mail, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))


# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.Fast)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    print("Liste paires: ", listPairs)

    # Some checks
    test.assertSequenceEqual(
        listPairs,
        [
            (0, 19),
            (0, 24),
            (3, 19),
            (3, 14),
            (3, 9),
            (1, 24),
            (1, 25),
            (1, 19),
            (1, 26),
            (1, 20),
            (1, 21),
            (6, 9),
            (4, 9),
            (4, 14),
            (4, 10),
            (4, 19),
            (4, 15),
            (4, 11),
            (4, 20),
            (4, 16),
            (4, 21),
            (2, 21),
            (2, 26),
            (2, 22),
            (2, 27),
            (7, 9),
            (7, 10),
            (7, 11),
            (5, 21),
            (5, 22),
            (5, 16),
            (5, 17),
            (5, 11),
            (5, 12),
            (8, 11),
            (8, 12),
        ],
    )

# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.Legacy)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
print("Nb pairs: ", nbPairs)

if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    print("Liste paires: ", listPairs)

    # Some checks
    test.assertSequenceEqual(
        listPairs,
        [
            (0, 19),
            (0, 24),
            (3, 19),
            (3, 14),
            (3, 9),
            (1, 24),
            (1, 25),
            (1, 19),
            (1, 26),
            (1, 20),
            (1, 21),
            (6, 9),
            (4, 9),
            (4, 14),
            (4, 10),
            (4, 19),
            (4, 15),
            (4, 11),
            (4, 20),
            (4, 16),
            (4, 21),
            (2, 21),
            (2, 26),
            (2, 22),
            (2, 27),
            (7, 9),
            (7, 10),
            (7, 11),
            (5, 21),
            (5, 22),
            (5, 16),
            (5, 17),
            (5, 11),
            (5, 12),
            (8, 11),
            (8, 12),
        ],
    )

# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.BrutForce)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    print("Liste paires: ", listPairs)

    test.assertSequenceEqual(
        listPairs,
        [
            (0, 19),
            (0, 24),
            (1, 19),
            (1, 20),
            (1, 21),
            (1, 24),
            (1, 25),
            (1, 26),
            (2, 21),
            (2, 22),
            (2, 26),
            (2, 27),
            (3, 9),
            (3, 14),
            (3, 19),
            (4, 9),
            (4, 10),
            (4, 11),
            (4, 14),
            (4, 15),
            (4, 16),
            (4, 19),
            (4, 20),
            (4, 21),
            (5, 11),
            (5, 12),
            (5, 16),
            (5, 17),
            (5, 21),
            (5, 22),
            (6, 9),
            (7, 9),
            (7, 10),
            (7, 11),
            (8, 11),
            (8, 12),
        ],
    )

# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(1)
meshPair.setPair("CONT_HAUT", "CONT_BAS")
# meshPair.setMethod("NEW")
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()


if nbPairs != 0:
    listPairs = meshPair.getListOfPairs()
    print("Liste paires: ", listPairs)

    # Some checks
    test.assertSequenceEqual(
        listPairs,
        [
            (9, 3),
            (9, 6),
            (9, 4),
            (9, 7),
            (14, 4),
            (14, 3),
            (10, 7),
            (10, 4),
            (19, 3),
            (19, 4),
            (19, 0),
            (19, 1),
            (15, 4),
            (11, 4),
            (11, 7),
            (11, 5),
            (11, 8),
            (24, 1),
            (24, 0),
            (20, 1),
            (20, 4),
            (16, 4),
            (16, 5),
            (12, 8),
            (12, 5),
            (25, 1),
            (21, 4),
            (21, 5),
            (21, 1),
            (21, 2),
            (17, 5),
            (26, 1),
            (26, 2),
            (22, 2),
            (22, 5),
            (27, 2),
        ],
    )

FIN()
