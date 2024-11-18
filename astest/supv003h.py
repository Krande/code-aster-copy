# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

import sys


from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import ExecutionParameter, Options

CA.init("--test", ERREUR=_F(ERREUR_F="EXCEPTION", ALARME="EXCEPTION"))
test = CA.TestCase()

MAIL = LIRE_MAILLAGE(FORMAT="MED")

MODELE = AFFE_MODELE(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

MAT = DEFI_MATERIAU(ELAS=_F(E=204000000000.0, NU=0.3, RHO=7800.0))

CHMAT = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MAT))

BLOCAGE = AFFE_CHAR_MECA(MODELE=MODELE, DDL_IMPO=_F(GROUP_MA="BASE", DX=0.0, DY=0.0, DZ=0.0))

size = 3
list_vect = []
for i in range(size):
    list_vect.append({"VECTEUR": CO(f"tmp{i}"), "OPTION": "CHAR_MECA", "CHARGE": BLOCAGE})


ASSEMBLAGE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    CHARGE=BLOCAGE,
    NUME_DDL=CO("NUMEDDL"),
    MATR_ASSE=(
        _F(MATRICE=CO("RIGIDITE"), OPTION="RIGI_MECA"),
        _F(MATRICE=CO("MASSE"), OPTION="MASS_MECA"),
    ),
    VECT_ASSE=list_vect,
)

# on vérifie que le premier et dernier vect_asse existent
test.assertTrue(tmp0.getType(), "VECT_ASSE")
last = globals()[f"tmp{size-1}"]
test.assertTrue(last.getType(), "VECT_ASSE")

ExecutionParameter().disable(Options.UseLegacyMode)
try:
    restup = ASSEMBLAGE(
        MODELE=MODELE,
        CHAM_MATER=CHMAT,
        CHARGE=BLOCAGE,
        NUME_DDL=CO("NUMEDDL"),
        MATR_ASSE=(
            _F(MATRICE=CO("RIGIDITE"), OPTION="RIGI_MECA"),
            _F(MATRICE=CO("MASSE"), OPTION="MASS_MECA"),
        ),
        VECT_ASSE=list_vect,
    )
except CA.AsterError as exc:
    test.assertEqual(exc.id_message, "SUPERVIS2_90")
    test.assertTrue(sys.version_info < (3, 7))
else:
    # main (None) + NUMEDDL + RIGIDITE + MASSE + tmp*
    test.assertEqual(len(restup), 1 + 3 + size)
    test.assertIsNone(restup.main, msg="restup.main should be None")
    test.assertIsNotNone(restup.NUMEDDL, msg="restup.NUMEDDL should exist")
    test.assertIsNotNone(restup.MASSE, msg="restup.MASSE should exist")
    test.assertIsNotNone(restup.RIGIDITE, msg="restup.RIGIDITE should exist")

test.printSummary()

CA.close()
