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


E_REF = 5.e11
E = DEFI_FONCTION(INFO=1,
                  INTERPOL='LIN',
                  NOM_PARA='TEMP',
                  NOM_RESU='TOUTRESU',
                  PROL_DROITE='CONSTANT',
                  PROL_GAUCHE='CONSTANT',
                  VALE=(0.0, E_REF, 1.0, E_REF),
                  VERIF='CROISSANT')


NU = DEFI_CONSTANTE(VALE=0.3)
RHO = DEFI_CONSTANTE(VALE=9800.)
ALPHA = DEFI_CONSTANTE(VALE=1.e-6)

acier = DEFI_MATERIAU(ELAS_FO=_F(E=E,
                               NU=NU,
                               RHO=RHO,
                               ALPHA=ALPHA,
                               TEMP_DEF_ALPHA=0.0,),)



# nombre de pas de temps
npas = 2
dt = 1.0 / npas
CHTEMP_VAR = []

PI = pi

FCHAR = FORMULE(VALE='cos(Z/PI)',
                NOM_PARA='Z',
                PI=PI)

CHGEOM = CREA_CHAMP(TYPE_CHAM='NOEU_GEOM_R',
                    OPERATION='EXTR',
                    MAILLAGE=mesh,
                    NOM_CHAM='GEOMETRIE',)

CHFONC = CREA_CHAMP(TYPE_CHAM='NOEU_NEUT_F',
                    OPERATION='AFFE',
                    MODELE=model,
                    AFFE=_F(TOUT='OUI',
                            NOM_CMP='X1',
                            VALE_F=FCHAR,),)

CHTEMPF = CREA_CHAMP(TYPE_CHAM='NOEU_NEUT_R',
                     OPERATION='EVAL',
                     CHAM_F=CHFONC,
                     CHAM_PARA=CHGEOM,)

CHTEMP = CREA_CHAMP(TYPE_CHAM='NOEU_TEMP_R',
                    OPERATION='ASSE',
                    MODELE=model,
                    ASSE=_F(TOUT='OUI',
                            CHAM_GD=CHTEMPF,
                            NOM_CMP='X1',
                            NOM_CMP_RESU='TEMP',),)

# on crée les différents champs de TEMP
for i in range(npas+1):
    time = i * dt
    new_temp = CHTEMP
    CHTEMP_VAR.append({'INST': time, 'CHAM_GD': new_temp})

TEMP = CREA_RESU(OPERATION='AFFE',
                 TYPE_RESU='EVOL_THER',
                 NOM_CHAM='TEMP',
                 AFFE=CHTEMP_VAR,
                 )


mater = AFFE_MATERIAU(MAILLAGE=mesh,
                      AFFE=_F(TOUT='OUI',
                              MATER=acier,),
                     AFFE_VARC=_F(TOUT='OUI',
                                   NOM_VARC='TEMP',
                                   EVOL=TEMP,
                                   VALE_REF=0.,),)


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
SOLUT = MECA_STATIQUE(MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      LIST_INST=LIST,
                      INFO=1)

SOLUN = MECA_STATIQUE2(MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      LIST_INST=LIST,
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


# =========================================================
#            REALISATION DES TESTS
# =========================================================

DIF_DEPL = DEPL_REF - DEPL

TEST_RESU(CHAM_NO=(_F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_DEPL,
                        VALE_CALC=1.5063505998114124E-12,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_DEPL,
                        VALE_CALC=1.7053025658242404E-12,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     ),
          )

FIN()
