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
import numpy as np

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

MA = LIRE_MAILLAGE(FORMAT="ASTER", UNITE=20)

MO = AFFE_MODELE(MAILLAGE=MA, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="COQUE_AXIS"))

BETON = DEFI_MATERIAU(ELAS=_F(E=3.7272000000e10, NU=0.0, RHO=2400.0))

COQUE = AFFE_CARA_ELEM(MODELE=MO, INFO=1, COQUE=_F(GROUP_MA=("COQUE"), EPAIS=0.5, COQUE_NCOU=4))

CHMAT = AFFE_MATERIAU(MAILLAGE=MA, AFFE=_F(TOUT="OUI", MATER=BETON))

BLOCAGE = AFFE_CHAR_MECA(MODELE=MO, DDL_IMPO=(_F(GROUP_NO="ENC", DX=0.0, DY=0.0, DRZ=0.0),))
CHARGE = AFFE_CHAR_MECA(MODELE=MO, FORCE_NODALE=(_F(GROUP_NO="CHA", FX=0, FY=-100)))

FOFO = DEFI_FONCTION(
    NOM_PARA="INST", VALE=(0.0, 0.0, 3.5, 3.5), PROL_DROITE="EXCLU", PROL_GAUCHE="EXCLU"
)


LINST = DEFI_LIST_REEL(
    DEBUT=0.0,
    INTERVALLE=(
        _F(JUSQU_A=0.1, NOMBRE=2),
        _F(JUSQU_A=1.4, NOMBRE=10),
        _F(JUSQU_A=2.84, NOMBRE=9),
        _F(JUSQU_A=3.0, NOMBRE=10),
    ),
)

U2 = MECA_STATIQUE(
    MODELE=MO,
    CHAM_MATER=CHMAT,
    CARA_ELEM=COQUE,
    EXCIT=(_F(CHARGE=BLOCAGE), _F(CHARGE=CHARGE, FONC_MULT=FOFO)),
    LIST_INST=LINST,
)

###  Tester les opérateurs arithmétques  et quelques méthodes ###
###                 de FieldsOnCells class                    ###


fieldOnNode1 = U2.getField("DEPL", 32)
fieldOnElem1 = U2.getField("SIEF_ELGA", 32)
fieldOnElem2 = U2.getField("SIEF_ELGA", 31)

# subscript operator __getitem__
valElem1 = fieldOnElem1.getValues()
test.assertAlmostEqual(fieldOnElem1[0], valElem1[0])

# __len__ and size()
test.assertAlmostEqual(len(fieldOnElem1), 192)
test.assertAlmostEqual(fieldOnElem1.size(), 192)


fieldOnElem3 = fieldOnElem2 - fieldOnElem1
test.assertAlmostEqual(fieldOnElem3[0], 1.742837573202337)

fieldOnElem3 = fieldOnElem2 + fieldOnElem1
test.assertAlmostEqual(fieldOnElem3[0], -651.8212523776659)

fieldOnElem3 = -fieldOnElem1
test.assertAlmostEqual(fieldOnElem3[0], -fieldOnElem1[0])

fieldOnElem3 = fieldOnElem1.copy()
fieldOnElem3 += fieldOnElem1
test.assertAlmostEqual(fieldOnElem3[0], 2 * fieldOnElem1[0])

fieldOnElem3 = fieldOnElem2.copy()
fieldOnElem3 -= fieldOnElem1

# this force "fieldOnElem2._values->updateValuePointer()"
valElem2 = fieldOnElem2.getValues()

test.assertAlmostEqual(fieldOnElem3[0], valElem2[0] - valElem1[0])

fieldOnElem3 = 10 * fieldOnElem2
test.assertAlmostEqual(fieldOnElem3[0], 10 * fieldOnElem2[0])

# __setitem__
fieldOnElem3[0] = 100
test.assertAlmostEqual(fieldOnElem3[0], 100)

# transform method

import math


def myfunc(x):
    return math.sin(x)


def myfuncMax(x):
    return min(1.0, max(0, x))


def myfuncErr(x):
    return str(x)


fieldOnElem4 = fieldOnElem3.transform(myfunc)
test.assertAlmostEqual(fieldOnElem4[0], math.sin(100))

fieldOnNode2 = fieldOnNode1.transform(myfuncMax)
test.assertAlmostEqual(fieldOnNode2[0], 0.0)
test.assertAlmostEqual(fieldOnNode2[1], 1.0)

with test.assertRaises(CA.AsterError):
    fieldOnNode3 = fieldOnNode1.transform(myfuncErr)
    fieldOnElem5 = fieldOnElem3.transform(myfuncErr)

# norms

vals = np.array(fieldOnElem1.getValues())
test.assertAlmostEqual(fieldOnElem1.norm("NORM_2"), np.sqrt(vals.dot(vals)))
test.assertAlmostEqual(fieldOnElem1.norm("NORM_1"), np.sum(np.abs(vals)))
test.assertAlmostEqual(fieldOnElem1.norm("NORM_INFINITY"), np.max(np.abs(vals)))

# dot
test.assertAlmostEqual(fieldOnElem1.dot(fieldOnElem1), vals.dot(vals))

test.printSummary()

FIN()
