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

DEBUT(CODE=_F(NIV_PUB_WEB='INTERNET'),)

test = code_aster.TestCase()

rank = code_aster.MPI.COMM_WORLD.Get_rank()

MA = code_aster.ParallelMesh()
MA.readMedFile("zzzz504y/%d.med" % rank, True)
# MA=LIRE_MAILLAGE()

MOD = AFFE_MODELE(MAILLAGE=MA,
                  AFFE=(_F(GROUP_MA=('Beton', 'Encast', 'Press',),
                           PHENOMENE='MECANIQUE',
                           MODELISATION='3D',),
                        _F(GROUP_MA='Cable0',
                           PHENOMENE='MECANIQUE',
                           MODELISATION='BARRE',),
                        ),)


# Orientation de tous les elements de surface
MA = MODI_MAILLAGE(reuse=MA,
                   MAILLAGE=MA,
                   ORIE_PEAU_3D=_F(GROUP_MA='Press'),
                   )

# Definition et effectation des materiaux

BTN_GEN = DEFI_MATERIAU(ELAS=_F(E=30000e6,
                                NU=0.2,
                                RHO=2500.,),
                        BPEL_BETON=_F(PERT_FLUA=0.0, PERT_RETR=0.0),)

ACI_CAB = DEFI_MATERIAU(ELAS=_F(E=200000e6,
                                NU=0.0,
                                RHO=7800,),
                        BPEL_ACIER=_F(RELAX_1000=0., F_PRG=0.),)


MATER = AFFE_MATERIAU(MAILLAGE=MA,
                      AFFE=(_F(GROUP_MA='Beton',
                               MATER=BTN_GEN,),
                            _F(GROUP_MA='Cable0',
                               MATER=ACI_CAB,),
                            ),)
# Definition et affectation des caracteristiques des elements de structure
ELEM = AFFE_CARA_ELEM(MODELE=MOD,
                      BARRE=_F(GROUP_MA='Cable0',
                               SECTION='CERCLE',
                               CARA='R',
                               VALE=0.005,),
                      )

LIAISON = AFFE_CHAR_MECA(MODELE=MOD,
                         DOUBLE_LAGRANGE='OUI',
                         LIAISON_MAIL=_F(GROUP_MA_MAIT='Beton', GROUP_MA_ESCL='Cable0')
                         )

GRAV = AFFE_CHAR_MECA(MODELE=MOD,
                      DOUBLE_LAGRANGE='OUI',
                      PESANTEUR=_F(GRAVITE=9.81, DIRECTION=(0., 0., -1.),),)

APP = AFFE_CHAR_CINE(MODELE=MOD,
                     MECA_IMPO=_F(GROUP_MA='Encast',
                                  DX=0.0,
                                  DY=0.0,
                                  DZ=0.0,),)

PRESS = AFFE_CHAR_MECA(MODELE=MOD,
                       DOUBLE_LAGRANGE='OUI',
                       PRES_REP=(_F(GROUP_MA='Press',
                                    PRES=1e6),
                                 ),
                       )
RAMPE = DEFI_FONCTION(NOM_PARA='INST',
                      VALE=(0., 0., 1.0, 1.0),
                      PROL_DROITE='CONSTANT',
                      )


INST = DEFI_LIST_REEL(DEBUT=0.,
                      INTERVALLE=_F(JUSQU_A=1,
                                    NOMBRE=1,),)

RESU = MECA_STATIQUE(MODELE=MOD,
                     CHAM_MATER=MATER,
                     CARA_ELEM=ELEM,
                     LIST_INST=INST,
                     EXCIT=(_F(CHARGE=APP),
                            _F(CHARGE=PRESS, FONC_MULT=RAMPE),
                            _F(CHARGE=LIAISON),),
                     OPTION='SANS',
                     INFO=2,)


TEST_RESU(
    RESU=_F(
        CRITERE='ABSOLU',
        GROUP_NO='Noeud1',
        NOM_CHAM='DEPL',
        NOM_CMP='DY',
        NUME_ORDRE=2,
        PRECISION=1.e-6,
        REFERENCE='AUTRE_ASTER',
        RESULTAT=RESU,
        VALE_CALC=-0.000505974,
        VALE_REFE=-0.000505974,
    ))


test.assertTrue(True)
test.printSummary()

FIN()
