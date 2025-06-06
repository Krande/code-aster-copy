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

test = CA.TestCase()

# extract from zzzz351a
DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))

MAIL = CA.Mesh()
MAIL.readMedFile("zzzz351a.mmed")

MODELE = AFFE_MODELE(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

MAT = DEFI_MATERIAU(ELAS=_F(E=200.0e9, NU=0.3, RHO=8000.0))

CHMAT = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MAT))

BLOCAGE = AFFE_CHAR_MECA(MODELE=MODELE, DDL_IMPO=_F(GROUP_MA="BLOK", DX=0.0, DY=0.0, DZ=0.0))

elem = CALC_MATR_ELEM(CHAM_MATER=CHMAT, CHARGE=BLOCAGE, MODELE=MODELE, OPTION="RIGI_MECA")

test.assertEqual(elem.getType(), "MATR_ELEM_DEPL_R")

FIN()

test.printSummary()
