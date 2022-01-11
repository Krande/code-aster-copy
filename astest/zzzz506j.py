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
from code_aster.MacroCommands.NonLinearSolver import NonLinearSolver, TimeStepper

DEBUT(CODE=_F(NIV_PUB_WEB='INTERNET',),
      DEBUG=_F(SDVERI='OUI',), INFO=1,)


mesh = LIRE_MAILLAGE(FORMAT='MED', UNITE=20)

model = AFFE_MODELE(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                            PHENOMENE='MECANIQUE',
                            MODELISATION='3D',),)

FTRACTUB = DEFI_FONCTION(NOM_PARA='EPSI',
                         VALE=(1000.,  2000.,
                               2000.,  5000.,),
                         PROL_DROITE='LINEAIRE',
                         PROL_GAUCHE='LINEAIRE',)

acier = DEFI_MATERIAU(ELAS=_F(E=200000.,
                              NU=0.3,),
                      ECRO_LINE=_F(D_SIGM_EPSI=2000.,
                                   SY=200.,),)

mater = AFFE_MATERIAU(MAILLAGE=mesh,
                      AFFE=_F(TOUT='OUI',
                              MATER=acier,),)


encast = AFFE_CHAR_MECA(MODELE=model,
                        DDL_IMPO=(_F(GROUP_MA='BAS',
                                      DX=0, DY=0.0, DZ=0.0,
                                      ),),)

depl = AFFE_CHAR_MECA(MODELE=model,
                      DDL_IMPO=(_F(GROUP_MA='HAUT',
                                    DZ=1.0,
                                    ),),)

LIST = DEFI_LIST_REEL(DEBUT=0.0,
                      INTERVALLE=_F(JUSQU_A=1.0,
                                    NOMBRE=2,),)

RAMPE = DEFI_FONCTION(NOM_PARA='INST',
                      VALE=(0., 0.,
                            1000., 1000.))

# STAT_NON_LINE DE REFERENCE
SOLUT = STAT_NON_LINE(MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      COMPORTEMENT=_F(RELATION='VMIS_ISOT_LINE',
                                      DEFORMATION="GDEF_LOG"),
                      NEWTON=_F(REAC_INCR=1,
                                PREDICTION='ELASTIQUE',
                                MATRICE='TANGENTE',
                                REAC_ITER=1,),
                      CONVERGENCE=_F(RESI_GLOB_RELA=1e-8,),
                      INCREMENT=_F(LIST_INST=LIST,),
                      INFO=1)

SOLUN = MECA_NON_LINE(MODELE=model,
                      CHAM_MATER=mater,
                      EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE),
                             _F(CHARGE=depl, FONC_MULT=RAMPE),),
                      COMPORTEMENT=_F(RELATION='VMIS_ISOT_LINE',
                                      DEFORMATION="GDEF_LOG"),
                      NEWTON=_F(REAC_INCR=1,
                                PREDICTION='ELASTIQUE',
                                MATRICE='TANGENTE',
                                REAC_ITER=1,),
                      CONVERGENCE=_F(RESI_GLOB_RELA=1e-8,),
                      INCREMENT=_F(LIST_INST=LIST,),
                      INFO=1)

# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant
SIGMA_REF = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='ELGA_SIEF_R',
                       NOM_CHAM='SIEF_ELGA', RESULTAT=SOLUT,
                       INST=1.
                       )

VARI_REF = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='ELGA_VARI_R',
                      NOM_CHAM='VARI_ELGA', RESULTAT=SOLUT,
                      INST=1.
                      )


SIGMA = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='ELGA_SIEF_R',
                   NOM_CHAM='SIEF_ELGA', RESULTAT=SOLUN,
                   INST=1.
                   )

VARI = CREA_CHAMP(OPERATION='EXTR', TYPE_CHAM='ELGA_VARI_R',
                  NOM_CHAM='VARI_ELGA', RESULTAT=SOLUN,
                  INST=1.
                  )

# =========================================================
#            REALISATION DES TESTS
# =========================================================

DIF_SIG = SIGMA_REF - SIGMA

DIF_VAR = VARI_REF - VARI

TEST_RESU(CHAM_ELEM=(_F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_SIG,
                        VALE_CALC=1.5063505998114124E-12,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        PRECISION=1.E-08,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_SIG,
                        VALE_CALC=1.7053025658242404E-12,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        ORDRE_GRANDEUR=5.E-03,
                        PRECISION=1.E-08,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_VAR,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='AUTRE_ASTER',
                        ORDRE_GRANDEUR=5.0e-3,
                        PRECISION=1.E-08,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_VAR,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI',),
                     ),
          )

FIN()
