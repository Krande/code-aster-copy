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


from code_aster.Commands import *
from code_aster import CA

CA.init("--test", "--continue")

test = CA.TestCase()

test.assertTrue("BLOK" in MAIL.getGroupsOfCells(), msg="BLOK in groups")
test.assertTrue("VOL" in MAIL.getGroupsOfCells(), msg="VOL in groups")

support = MODELE.getMesh()
test.assertTrue("VOL" in support.getGroupsOfCells(), msg="VOL in groups from model")
del support

asse = ASSEMBLAGE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    CHARGE=BLOCAGE,
    NUME_DDL=CO("NUMEDDL"),
    MATR_ASSE=(_F(MATRICE=CO("K1"), OPTION="RIGI_MECA"), _F(MATRICE=CO("M1"), OPTION="MASS_MECA")),
)
test.assertIsNone(asse, msg="returns nothing")

# check for dependencies
deps = K1.getDependencies()
test.assertEqual(len(deps), 0, msg="K1 dependencies")

# 1. Calcul de reference avec les matrices "completes" :
# --------------------------------------------------------
MODE1 = CALC_MODES(
    OPTION="BANDE",
    MATR_RIGI=K1,
    MATR_MASS=M1,
    CALC_FREQ=_F(FREQ=(10.0, 250.0)),
    VERI_MODE=_F(SEUIL=1e-03),
)

TEST_RESU(RESU=_F(RESULTAT=MODE1, NUME_MODE=2, PARA="FREQ", VALE_CALC=85.631015163879))

test.printSummary()

FIN()
