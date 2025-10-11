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
meshPair.setVerbosity(2)
meshPair.setPair("CONT_BAS", "CONT_HAUT")
meshPair.setMethod(PairingMethod.Fast)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
listPairs = meshPair.getListOfPairs()

# Get quadrature points
indexPair = 4
quadPoin = meshPair.getQuadraturePoints(indexPair)
print(" => Coordinates of quadrature points (global space):", quadPoin)

# Some checks
test.assertSequenceEqual(
    listPairs,
    [
        (465, 314),
        (465, 316),
        (465, 332),
        (465, 331),
        (465, 330),
        (468, 330),
        (468, 313),
        (468, 331),
        (468, 328),
        (468, 312),
        (466, 330),
        (466, 331),
        (466, 360),
        (466, 332),
        (466, 362),
        (466, 314),
        (466, 344),
        (466, 364),
        (466, 316),
        (466, 345),
        (466, 358),
        (466, 335),
        (466, 357),
        (471, 312),
        (469, 312),
        (469, 328),
        (469, 313),
        (469, 352),
        (469, 330),
        (469, 341),
        (469, 356),
        (469, 360),
        (469, 343),
        (469, 347),
        (469, 362),
        (469, 364),
        (467, 357),
        (467, 359),
        (467, 358),
        (467, 355),
        (467, 365),
        (467, 364),
        (467, 354),
        (467, 363),
        (467, 351),
        (467, 361),
        (467, 353),
        (472, 312),
        (472, 328),
        (472, 352),
        (472, 341),
        (472, 343),
        (472, 327),
        (470, 364),
        (470, 362),
        (470, 365),
        (470, 360),
        (470, 363),
        (470, 356),
        (470, 350),
        (470, 361),
        (470, 347),
        (470, 348),
        (470, 353),
        (470, 343),
        (470, 333),
        (470, 349),
        (470, 327),
        (470, 326),
        (473, 327),
        (473, 343),
        (473, 333),
        (473, 326),
    ],
)


# Generate pairs
meshPair = CA.MeshPairing()
meshPair.setMesh(Mail)
meshPair.setVerbosity(2)
meshPair.setPair("CONT_HAUT", "CONT_BAS")
meshPair.setMethod(PairingMethod.Fast)
meshPair.compute()

# Get pairs
nbPairs = meshPair.getNumberOfPairs()
listPairs = meshPair.getListOfPairs()

# Some checks
test.assertSequenceEqual(
    listPairs,
    [
        (312, 468),
        (312, 471),
        (312, 469),
        (312, 472),
        (328, 472),
        (328, 469),
        (328, 468),
        (313, 468),
        (313, 469),
        (352, 469),
        (352, 472),
        (330, 469),
        (330, 468),
        (330, 466),
        (330, 465),
        (341, 472),
        (341, 469),
        (356, 469),
        (356, 470),
        (331, 465),
        (331, 468),
        (331, 466),
        (360, 466),
        (360, 469),
        (360, 470),
        (343, 469),
        (343, 472),
        (343, 470),
        (343, 473),
        (347, 470),
        (347, 469),
        (332, 466),
        (332, 465),
        (362, 470),
        (362, 469),
        (362, 466),
        (327, 473),
        (327, 472),
        (327, 470),
        (350, 470),
        (314, 465),
        (314, 466),
        (344, 466),
        (364, 466),
        (364, 469),
        (364, 467),
        (364, 470),
        (333, 470),
        (333, 473),
        (363, 470),
        (363, 467),
        (348, 470),
        (316, 466),
        (316, 465),
        (345, 466),
        (358, 467),
        (358, 466),
        (365, 470),
        (365, 467),
        (326, 473),
        (326, 470),
        (361, 467),
        (361, 470),
        (349, 470),
        (335, 466),
        (357, 466),
        (357, 467),
        (359, 467),
        (351, 467),
        (353, 470),
        (353, 467),
        (355, 467),
        (354, 467),
    ],
)


FIN()
