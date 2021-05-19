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

import numpy as N
import code_aster
from code_aster.Commands import *
from code_aster import MPI

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.COMM_WORLD.Get_rank()

MAIL= LIRE_MAILLAGE(FORMAT='MED',
                       PARTITIONNEUR='PTSCOTCH',
                       INFO=1,
                       )
DEFI_GROUP(reuse=MAIL, MAILLAGE=MAIL, CREA_GROUP_NO=_F(TOUT_GROUP_MA='OUI',))

MATER=DEFI_MATERIAU(ELAS=_F(E=10000.0,
                            NU=0.,
                            RHO=1.0,),
                            );

affectMat = code_aster.MaterialField(MAIL)
affectMat.addMaterialsOnMesh(MATER)
affectMat.buildWithoutExternalVariable()

MODT=AFFE_MODELE(MAILLAGE=MAIL,
                 AFFE=_F(GROUP_MA=('S11',    'S31', 'S12',     'S32'),
                         PHENOMENE='MECANIQUE',
                         MODELISATION='D_PLAN',),)


CHT1 = AFFE_CHAR_MECA(MODELE=MODT,
                      PESANTEUR=_F(GRAVITE=1.0,
                                   DIRECTION=(0.0, -1.0, 0.0),),
                      INFO=1,
                      VERI_NORM='NON',)

charCine = AFFE_CHAR_CINE(MODELE=MODT,
            MECA_IMPO=(_F(GROUP_MA=("Bas1"), DX=0., DY=1.),
                       _F(GROUP_MA=("Bas3"), DX=0., DY=-1.),))


resu=MECA_STATIQUE( MODELE=MODT,
                        CHAM_MATER=affectMat,
                        EXCIT  =(_F( CHARGE = CHT1), _F( CHARGE = charCine)),
                       )


# resu.printMedFile('/tmp/test2_{}.resu.med'.format(rank))

TEST_RESU(
    RESU=_F(
        CRITERE='ABSOLU',
        GROUP_NO='Point1',
        NOM_CHAM='DEPL',
        NOM_CMP='DY',
        NUME_ORDRE=1,
        PRECISION=1.e-6,
        REFERENCE='AUTRE_ASTER',
        RESULTAT=resu,
        VALE_CALC=0.999549999999998,
        VALE_REFE=0.999549999999998,
    ))

TEST_RESU(
    RESU=_F(
        CRITERE='ABSOLU',
        GROUP_NO='Point4',
        NOM_CHAM='DEPL',
        NOM_CMP='DY',
        NUME_ORDRE=1,
        PRECISION=1.e-6,
        REFERENCE='AUTRE_ASTER',
        RESULTAT=resu,
        VALE_CALC=-1.0004499999999994,
        VALE_REFE=-1.0004499999999994,
    ))


test.printSummary()

FIN()
