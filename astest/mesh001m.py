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

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION", ERREUR_F="EXCEPTION"))

test = CA.TestCase()

fmt_raison = (
    "-" * 80
    + """

   Exception interceptee
   Message : %s

"""
    + "-" * 80
    + "\n"
)

is_ok = 0
try:
    mesh_3d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=20)
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "MODELISA4_10":
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

is_ok = 0
try:
    mesh_2d = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=21)
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "MODELISA4_10":
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

test.printSummary()

FIN()
