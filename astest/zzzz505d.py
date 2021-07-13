# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

DEBUT(CODE=_F(NIV_PUB_WEB='INTERNET',),
      DEBUG=_F(SDVERI='OUI',),
      INFO=1,)

test = code_aster.TestCase()

mesh=LIRE_MAILLAGE(FORMAT='MED', UNITE=20)

model=AFFE_MODELE(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                         PHENOMENE='MECANIQUE',
                         MODELISATION='3D',),)
# Very high elasticity limit to simulate elasticity
acier=DEFI_MATERIAU(ELAS=_F(E=200000.,
                            NU=0.3,),
                    ECRO_LINE=_F(D_SIGM_EPSI=2000.,
                                 SY=200000.,),)

mater=AFFE_MATERIAU(MAILLAGE=mesh,
                   AFFE=_F(TOUT='OUI',
                           MATER=acier,),)

# Build reference field with CREA_CHAMP

value=2.0

refe1=CREA_CHAMP(TYPE_CHAM='ELGA_SIEF_R',
           OPERATION='AFFE',
           MODELE=model,
           AFFE=_F(TOUT='OUI',
           NOM_CMP=('SIXX','SIYY','SIZZ',
                    'SIXY','SIYZ','SIXZ',),
           VALE=(value,value,value,value,value,value,),
                              ),)
refe2=CREA_CHAMP(OPERATION='AFFE',
                 TYPE_CHAM='ELGA_VARI_R',
                 MODELE=model,
                 PROL_ZERO='OUI',
                 AFFE=_F(TOUT='OUI', NOM_CMP= ('V1','V2'), VALE = (value,value),)
                 )

# Using Python binding for behaviour
study = code_aster.StudyDescription(model, mater)
dProblem = code_aster.DiscreteProblem(study)

# With default values: no initial state, no implex and info=1
behav = dProblem.createBehaviour(
    COMPORTEMENT=(
        _F(RELATION="VMIS_ISOT_LINE",
           TOUT="OUI"),))
# Build testfield with model and compor
testfield1=code_aster.FieldOnCellsReal(model, behav, "ELGA_SIEF_R")
testfield1.setValues(value)

testfield2=code_aster.FieldOnCellsReal(model, behav, "ELGA_VARI_R")
testfield2.setValues(value)

# Test
test.assertAlmostEqual(len(refe1.getValues()), len(testfield1.getValues()))
test.assertAlmostEqual(len(refe2.getValues()), len(testfield2.getValues()))
test.assertAlmostEqual(refe1.getValues(), testfield1.getValues())
test.assertAlmostEqual(refe2.getValues(), testfield2.getValues())
test.printSummary()

code_aster.close()
