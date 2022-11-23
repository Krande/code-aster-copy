# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

from code_aster import AsterError
from code_aster.Commands import *

DEBUT(CODE=_F(NIV_PUB_WEB="INTERNET"), DEBUG=_F(SDVERI="OUI"), ERREUR=_F(ERREUR_F="EXCEPTION"))

fmt_raison = (
    "-" * 80
    + """

   Exception interceptee
   Message : %s

"""
    + "-" * 80
    + "\n"
)


mesh0 = LIRE_MAILLAGE(FORMAT="MED")

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

is_ok = 0
try:
    model = AFFE_MODELE(
        AFFE=(
            _F(MODELISATION=("3D",), PHENOMENE="MECANIQUE", TOUT="OUI"),
            _F(
                MODELISATION=("3D_HHO",),
                FORMULATION="QUADRATIQUE",
                PHENOMENE="MECANIQUE",
                TOUT="OUI",
            ),
        ),
        MAILLAGE=mesh,
    )
except AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "MODELE1_10":
        is_ok = 1

TAB1 = CREA_TABLE(
    LISTE=(_F(PARA="TEST", TYPE_K="K8", LISTE_K="VALEUR  "), _F(PARA="BOOLEEN", LISTE_I=is_ok))
)
TEST_TABLE(
    REFERENCE="ANALYTIQUE",
    VALE_CALC_I=1,
    VALE_REFE_I=1,
    NOM_PARA="BOOLEEN",
    TABLE=TAB1,
    FILTRE=_F(NOM_PARA="TEST", VALE_K="VALEUR  "),
)

DETRUIRE(NOM=(TAB1))


FIN()
