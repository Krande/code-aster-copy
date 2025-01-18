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

from math import sqrt
from code_aster.Commands import *
from code_aster import CA

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

U2 = CALC_CHAMP(reuse=U2, RESULTAT=U2, CRITERES="SIEQ_ELGA")

last = U2.getLastIndex()
sigm = U2.getField("SIEF_ELGA", last).toSimpleFieldOnCells()
sieq = U2.getField("SIEQ_ELGA", last).toSimpleFieldOnCells()

test.assertAlmostEqual(sigm.getValue(0, 0, 0, 0), -326.7820449754335, 8)
test.assertIn(sigm.getPhysicalQuantity(), ("SIEF_R",))
test.assertIn(sigm.getLocalization(), ("ELGA",))
test.assertSequenceEqual(sigm.getComponents(), ["SIXX", "SIYY", "SIZZ", "SIXZ"])
test.assertEqual(sigm.getMaxNumberOfPoints(), 4)
test.assertEqual(sigm.getNumberOfPointsOfCell(0), 4)
test.assertEqual(sigm.getNumberOfCells(), 1)
test.assertEqual(sigm.getNumberOfComponents(), 4)
test.assertEqual(sigm.getNumberOfSubPointsOfCell(0), 12)

s1 = sieq.PRIN_1
s2 = sieq.PRIN_2
s3 = sieq.PRIN_3
vmis = sieq.VMIS
calc = (1.0 / sqrt(2.0)) * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s1 - s3) ** 2) ** 0.5

test.assertAlmostEqual(abs(vmis - calc).max(), 0.0, 8)

chcoor = CALC_CHAM_ELEM(MODELE=MO, CARA_ELEM=COQUE, OPTION="COOR_ELGA")
coor = chcoor.toSimpleFieldOnCells()

weight = coor.W
weight_sigm = weight.onSupportOf(s1)
# \int { PRIN_1 }
s1xw_c = s1 * weight_sigm
print(repr(s1xw_c))
test.assertAlmostEqual(s1xw_c.sum(), -2482.24915, 4)

test.printSummary()

FIN()
