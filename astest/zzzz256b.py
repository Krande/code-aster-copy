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

# Cas test qui vérifie la résolution d'un système linéaire non-symétrique
# voir issue31773

from code_aster.Commands import *
from code_aster import CA

test = CA.TestCase()

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))

MAIL = LIRE_MAILLAGE(FORMAT="MED", UNITE=20, PARTITIONNEUR="PTSCOTCH")

MAIL = MODI_MAILLAGE(
    reuse=MAIL,
    MAILLAGE=MAIL,
    ORIE_PEAU=_F(GROUP_MA_PEAU="IFS_2D", GROUP_MA_INTERNE=("s1Eau", "Eau")),
)

# -- CONSTRUCTION DU MODELE
MODELE = AFFE_MODELE(
    AFFE=(
        _F(GROUP_MA=("s1Tuyau", "Tuyau"), MODELISATION="3D", PHENOMENE="MECANIQUE"),
        _F(
            FORMULATION="U_P",
            GROUP_MA=("s1Eau", "Eau"),
            MODELISATION="3D_FLUIDE",
            PHENOMENE="MECANIQUE",
        ),
        _F(
            FORMULATION="U_P", GROUP_MA=("IFS_2D",), MODELISATION="FLUI_STRU", PHENOMENE="MECANIQUE"
        ),
    ),
    MAILLAGE=MAIL,
)


# -- DEFINITION DU MATERIAU
scaling = 1.0e10
acier = DEFI_MATERIAU(ELAS=_F(E=190.0e9 / scaling, NU=0.3, RHO=7800.0 / scaling))

eau = DEFI_MATERIAU(FLUIDE=_F(CELE_R=970.0, RHO=732.0 / scaling))

# -- AFFECTATION DES MATERIAUX
CHMAT = AFFE_MATERIAU(
    AFFE=(_F(GROUP_MA=("Tuyau"), MATER=(acier,)), _F(GROUP_MA=("IFS_2D", "Eau"), MATER=(eau,))),
    MODELE=MODELE,
)


# -- DEFINITION DES CONDITIONS AUX LIMITES
BLOQ = AFFE_CHAR_MECA(
    FACE_IMPO=(
        _F(DX=0.0, DY=0.0, DZ=0.0, GROUP_MA=("s1Tuyau",)),
        _F(GROUP_MA=("s1Eau",), PRES=0.0),
        _F(DX=0.0, DY=0.0, DZ=0.0, PRES=0.0, GROUP_MA=("Interface",)),
    ),
    MODELE=MODELE,
)

PESA = AFFE_CHAR_MECA(PESANTEUR=_F(GRAVITE=9.81, DIRECTION=(0.0, 0.0, -1.0)), MODELE=MODELE)


resu = MECA_STATIQUE(MODELE=MODELE, CHAM_MATER=CHMAT, EXCIT=(_F(CHARGE=BLOQ), _F(CHARGE=PESA)))

depl = resu.getField("DEPL", 1)
test.assertAlmostEqual(depl.norm("NORM_2"), 0.001181792186987481)


FIN()
