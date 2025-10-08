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

import time


# -------------------------------------
#       - Main
# -------------------------------------
DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

Mail = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

Mail = MODI_MAILLAGE(
    reuse=Mail, MAILLAGE=Mail, ORIE_PEAU=_F(GROUP_MA_PEAU=("CONT_HAUT", "CONT_BAS"))
)

MAT = DEFI_MATERIAU(ELAS=_F(E=20000, NU=0.3, ALPHA=0.01))

CHMAT = AFFE_MATERIAU(MAILLAGE=Mail, AFFE=_F(TOUT="OUI", MATER=MAT))

MODI = AFFE_MODELE(MAILLAGE=Mail, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"))

# -------------------------------------
#       - Slave side - CONT_BAS
# -------------------------------------
## -- OLD PAIRING METHOD
# - Generate pairs
meshPair_old = CA.MeshPairing()
meshPair_old.setMesh(Mail)
meshPair_old.setVerbosity(2)
meshPair_old.setPair("CONT_BAS", "CONT_HAUT")
meshPair_old.setMethod(PairingMethod.Legacy)

t_0_old = time.time()
meshPair_old.compute()
dur_old = time.time() - t_0_old

# - Get pairs
nbPairs_old = meshPair_old.getNumberOfPairs()
listPairs_old = meshPair_old.getListOfPairs()

# -- FAST PAIRING METHOD
# - Generate pairs
meshPair_fast = CA.MeshPairing()
meshPair_fast.setMesh(Mail)
meshPair_fast.setVerbosity(2)
meshPair_fast.setPair("CONT_BAS", "CONT_HAUT")
meshPair_fast.setMethod(PairingMethod.BrutForce)

t_0_fast = time.time()
meshPair_fast.compute()
dur_fast = time.time() - t_0_fast

# - Get pairs
nbPairs_fast = meshPair_fast.getNumberOfPairs()
listPairs_fast = meshPair_fast.getListOfPairs()

test.assertEqual(nbPairs_fast, nbPairs_old)


FIN()
