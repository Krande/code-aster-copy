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

import code_aster
from code_aster.Commands import *

code_aster.init("--test")

test = code_aster.TestCase()

MA=LIRE_MAILLAGE(FORMAT="ASTER",UNITE=20,)

MO=AFFE_MODELE( MAILLAGE=MA,
                AFFE=_F(TOUT='OUI', PHENOMENE='MECANIQUE', MODELISATION='COQUE_AXIS',),)

BETON=DEFI_MATERIAU(ELAS=_F(E=3.7272000000E10, NU=0.0,  RHO=2400.0,),)

COQUE=AFFE_CARA_ELEM(  MODELE=MO,
                        INFO=1,
                        COQUE = _F( GROUP_MA=('COQUE'),
                                    EPAIS = 0.5,
                                    COQUE_NCOU = 4,),
                    )

CHMAT=AFFE_MATERIAU(MAILLAGE=MA,
                    AFFE=_F(TOUT='OUI', MATER=BETON,),
                    )

BLOCAGE=AFFE_CHAR_MECA( MODELE=MO,
                        DDL_IMPO=(_F(GROUP_NO='ENC', DX=0.0, DY=0.0, DRZ=0.0,),
                                 ),)
CHARGE=AFFE_CHAR_MECA(  MODELE=MO,
                        FORCE_NODALE=(_F(GROUP_NO='CHA',FX = 0,FY = -100)
                                 ),)

FOFO=DEFI_FONCTION(NOM_PARA='INST',
                   VALE=(0.0,0.0,3.5,3.5),
                   PROL_DROITE='EXCLU',
                   PROL_GAUCHE='EXCLU',)


LINST=DEFI_LIST_REEL(   DEBUT=0.0,
                        INTERVALLE=(_F(JUSQU_A=0.1,  NOMBRE=2,),
                                    _F(JUSQU_A=1.4,  NOMBRE=10,),
                                    _F(JUSQU_A=2.84 ,  NOMBRE=9,),
                                    _F(JUSQU_A=3. ,  NOMBRE=10,),
                                    ),
                    )

U2 = MECA_STATIQUE(MODELE=MO,
                 CHAM_MATER=CHMAT,
                 CARA_ELEM=COQUE,
                 EXCIT=(_F(CHARGE=BLOCAGE,),
                        _F(CHARGE=CHARGE,
                           FONC_MULT=FOFO,),),
                 LIST_INST=LINST,
                 )

fieldOnElem = U2.getFieldOnCellsReal("SIEF_ELGA", 31)

sfon = fieldOnElem.exportToSimpleFieldOnCells()

test.assertAlmostEqual(sfon.getValue(0, 0, 0, 0), -325.03920740223253)
test.assertIn(sfon.getPhysicalQuantity(), ('SIEF_R',))
test.assertIn(sfon.getFieldLocation(), ('ELGA',))
test.assertSequenceEqual(sfon.getNameOfComponents(), ['SIXX', 'SIYY', 'SIZZ', 'SIXZ'])
test.assertEqual(sfon.getMaxNumberOfPoints(), 4)
test.assertEqual(sfon.getNumberOfPointsOfCell(0), 4)
test.assertEqual(sfon.getNumberOfCells(), 1)
test.assertEqual(sfon.getNumberOfComponents(), 4)
test.assertEqual(sfon.getNumberOfSubPointsOfCell(0), 12)

test.printSummary()

FIN()
