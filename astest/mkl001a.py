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


DEBUT(CODE=_F(NIV_PUB_WEB='INTERNET',),
      DEBUG=_F(SDVERI='OUI',), INFO=1,)

mesh = LIRE_MAILLAGE(FORMAT='MED', UNITE=20)

model = AFFE_MODELE(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                            PHENOMENE='MECANIQUE',
                            MODELISATION='3D',),)


acier = DEFI_MATERIAU(ELAS=_F(E=200000.,
                              NU=0.3,),)

mater = AFFE_MATERIAU(MAILLAGE=mesh,
                      AFFE=_F(TOUT='OUI',
                              MATER=acier,),)


encast = AFFE_CHAR_CINE(MODELE=model,
                        MECA_IMPO=(_F(GROUP_MA='BAS',
                                      DX=0, DY=0.0, DZ=0.0,
                                      ),),)

depl = AFFE_CHAR_CINE(MODELE=model,
                      MECA_IMPO=(_F(GROUP_MA='HAUT',
                                    DZ=1.0,
                                    ),),)

LIST = DEFI_LIST_REEL(DEBUT=0.0,
                      INTERVALLE=_F(JUSQU_A=1.0,
                                    NOMBRE=2,),)

RAMPE = DEFI_FONCTION(NOM_PARA='INST',
                      VALE=(0., 0.,
                            1000., 1000.))

# MECA_STATIQUE DE REFERENCE
SOLUT = MECA_STATIQUE2(MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      LIST_INST=LIST,
                      INFO=1)

print(SOLUT.getAccessParameters())

SOLUN = MECA_STATIQUE(MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      INST=0.0,
                      INFO=1)

print(SOLUN.getAccessParameters())

LIST2 = DEFI_LIST_REEL(VALE=(0.5, 1.0),)

SOLUN = MECA_STATIQUE(reuse=SOLUN, RESULTAT=SOLUN,
                      MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      LIST_INST=LIST2,
                      INFO=1)

# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant
DEPL_REF = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='NOEU_DEPL_R',
                       NOM_CHAM='DEPL', RESULTAT=SOLUT,
                       INST=1.
                       )

DEPL = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='NOEU_DEPL_R',
                   NOM_CHAM='DEPL', RESULTAT=SOLUN,
                   INST=1.
                   )

SIEF_REF = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='ELGA_SIEF_R',
                       NOM_CHAM='SIEF_ELGA', RESULTAT=SOLUT,
                       INST=1.
                       )

SIEF = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='ELGA_SIEF_R',
                       NOM_CHAM='SIEF_ELGA', RESULTAT=SOLUN,
                   INST=1.
                   )


# =========================================================
#            REALISATION DES TESTS
# =========================================================

DIF_DEPL = DEPL_REF - DEPL

DIF_SIEF = SIEF_REF - SIEF

TEST_RESU(CHAM_NO=(_F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        ORDRE_GRANDEUR=1e-8,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_DEPL,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        ORDRE_GRANDEUR=1e-8,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_DEPL,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     ),
          )

TEST_RESU(CHAM_ELEM=(_F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        ORDRE_GRANDEUR=1e-8,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_SIEF,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        ORDRE_GRANDEUR=1e-8,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_SIEF,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     ),
          )

FIN()
