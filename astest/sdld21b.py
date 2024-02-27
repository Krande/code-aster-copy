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

NU2 = NUME_DDL(MODELE=MO, CHARGE=(CL))

dt = 0.0081551968

tf = 255.0 * dt

LINST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=tf, PAS=dt))

CHA_ON = CREA_RESU(
    OPERATION="CONV_RESU",
    TYPE_RESU="EVOL_CHAR",
    CONV_RESU=_F(
        NUME_DDL=NU2,
        NOM_CHAM_INIT="ACCE",
        COEF=1.0e13,
        PRECISION=1.0e-6,
        CRITERE="RELATIF",
        LIST_INST=LINST,
        RESU_INIT=DNTPSAP2,
    ),
)

test.assertEqual(CHA_ON.getNumberOfIndexes(), 256)
test.assertAlmostEqual(CHA_ON.getTime(1), 0.0)
test.assertAlmostEqual(CHA_ON.getTime(2), dt)

test.printSummary()

FIN()
