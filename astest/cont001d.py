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
from code_aster.MacroCommands.defi_cont import DEFI_CONT

DEBUT(CODE="OUI", ERREUR=_F(ALARME="ALARME", ERREUR_F="EXCEPTION"), INFO=1)

test = CA.TestCase()

fmt_raison = (
    "-" * 80
    + """

   Exception interceptee
   Raison : %s

"""
    + "-" * 80
    + "\n"
)


Mail = LIRE_MAILLAGE(UNITE=20, PARTITIONNEUR="PTSCOTCH", FORMAT="MED")

MODI = AFFE_MODELE(MAILLAGE=Mail, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="D_PLAN"))

# Check no shared nodes
is_ok = 0
try:
    DEFICO = DEFI_CONT(MODELE=MODI, ZONE=(_F(GROUP_MA_MAIT="Group_2", GROUP_MA_ESCL="Group_2"),))
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "CONTACT1_1":
        is_ok = 1

test.assertEqual(is_ok, 1)

# check normals orientations
is_ok = 0
try:
    DEFIC1 = DEFI_CONT(MODELE=MODI, ZONE=(_F(GROUP_MA_MAIT="Group_1", GROUP_MA_ESCL="Group_2"),))
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "MODELISA4_24":
        is_ok = 1
test.assertEqual(is_ok, 1)


is_ok = 0
try:
    DEFIC2 = DEFI_CONT(
        MODELE=MODI, ZONE=(_F(GROUP_MA_MAIT="Group_1", VERI_NORM="NON", GROUP_MA_ESCL="Group_2"),)
    )
    is_ok = 1
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "MODELISA4_24":
        is_ok = 0

test.assertEqual(is_ok, 1)


# check mechanical model
is_ok = 0
try:
    MODI_THER = AFFE_MODELE(
        MAILLAGE=Mail, AFFE=_F(TOUT="OUI", PHENOMENE="THERMIQUE", MODELISATION="PLAN")
    )
    DEFIC3 = DEFI_CONT(
        MODELE=MODI_THER, ZONE=(_F(GROUP_MA_MAIT="Group_1", GROUP_MA_ESCL="Group_2"),)
    )
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "CONTACT1_2":
        is_ok = 1
test.assertEqual(is_ok, 1)

is_ok = 0
try:
    DEFICO = DEFI_CONT(
        MODELE=MODI,
        ZONE=(
            _F(GROUP_MA_MAIT="Group_1", VERI_NORM="NON", GROUP_MA_ESCL="Group_2"),
            _F(GROUP_MA_MAIT="Group_3", VERI_NORM="NON", GROUP_MA_ESCL="Group_2"),
        ),
    )
except CA.AsterError as err:
    print(fmt_raison % str(err))
    # on verifie que l'erreur fatale est bien celle que l'on attendait :
    if err.id_message == "CONTACT1_1":
        is_ok = 1

test.assertEqual(is_ok, 1)

FIN()
