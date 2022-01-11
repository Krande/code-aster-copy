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

# ***********************************************************************
#
# TITRE: ESSAI DE TRACTION AVEC LA LOI DE RANKINE
# MODIFICATION DU TEST SSNV515C
#
# ***********************************************************************

# ======================================================================
DEBUT(CODE=_F(NIV_PUB_WEB='INTERNET', ), );

# ***********************************************************************
#
#    MAILLAGE ET MODELE
#
# ***********************************************************************

MAILLAGE = LIRE_MAILLAGE(FORMAT='MED', PARTITIONNEUR='PTSCOTCH',);

MODELE = AFFE_MODELE(MAILLAGE=MAILLAGE,
                     AFFE=_F(TOUT='OUI',
                             PHENOMENE='MECANIQUE',
                             MODELISATION='AXIS', ), );

MAILLAGE = MODI_MAILLAGE(reuse=MAILLAGE,
                         MAILLAGE=MAILLAGE,
                         ORIE_PEAU=_F(GROUP_MA_PEAU=('DROIT', 'GAUCHE',
                                                     'BAS', 'HAUT'), ),
                         INFO=1, );

# ***********************************************************************
#
#    LISTE D'INSTANTS
#
# ***********************************************************************

tarret = 10.
tfin = 20.
npas = 30
temps_max = 30.

dtemps = temps_max / npas

ltemps = [dtemps * i for i in range(npas + 1)]

TEMPS = DEFI_LIST_REEL(DEBUT=0.,
                       INTERVALLE=(_F(JUSQU_A=20, NOMBRE=10, ),), );
TEMPS2 = DEFI_LIST_REEL(DEBUT=20,
                        INTERVALLE=(_F(JUSQU_A=temps_max, NOMBRE=5, ),), );

# ***********************************************************************
#
#    DONNEES MATERIAU
#
# ***********************************************************************

# modules mecaniques [kPa]
# ------------------------
# K=516.2E6
# G=238.2E6
# # =>
# YOUNG = 9.*K*G /(3.*K+G)
# POISSON = (3.*K-2.*G) /(6.*K+2.*G)

YOUNG = 1e+6
POISSON = 0.25
SIGMA_T = 1.E+3
K = YOUNG / 3. / (1. - 2. * POISSON)
G = YOUNG / 2. / (1. + POISSON)

SOL = DEFI_MATERIAU(ELAS=_F(E=YOUNG, NU=POISSON, ALPHA=0., ),
                    RANKINE=_F(SIGMA_T=SIGMA_T, ),
                    INFO=1, );

CHMAT = AFFE_MATERIAU(MAILLAGE=MAILLAGE,
                      AFFE=_F(TOUT='OUI', MATER=SOL, ), );

# ***********************************************************************
#
#    CHARGEMENTS
#
# ***********************************************************************

# pression de preconsolidation [en kPa]
P0 = 5.E+3
EPZZ = .03

npas = 300
dtemps = temps_max / npas
linst = [dtemps * i for i in range(npas)]

SIGLAT = AFFE_CHAR_MECA(MODELE=MODELE,
                        PRES_REP=_F(GROUP_MA=('DROIT',),
                                    PRES=P0, ), );

DEPHAUT = AFFE_CHAR_CINE(MODELE=MODELE,
                         MECA_IMPO=(_F(GROUP_MA=('HAUT',), DY=1., ),), );

DEPL_1 = AFFE_CHAR_CINE(MODELE=MODELE,
                        MECA_IMPO=(_F(GROUP_MA='BAS', DY=0., ),_F(GROUP_NO='B', DX=0., ),), );



COEF2 = DEFI_FONCTION(NOM_PARA='INST',
                      PROL_DROITE='CONSTANT',
                      VALE=(0., 0., 20, EPZZ,
                            temps_max, 0,), );

COEF3 = DEFI_FONCTION(NOM_PARA='INST',
                      PROL_DROITE='CONSTANT',
                      VALE=(0., 1.,), );


# ***********************************************************************
#
#    PRECONSOLIDATION INITIALE A 5KPA
#
# ***********************************************************************

SIG0 = CREA_CHAMP(INFO=2,
                  TYPE_CHAM='ELGA_SIEF_R',
                  OPERATION='AFFE',
                  MODELE=MODELE,
                  PROL_ZERO='OUI',
                  AFFE=_F(GROUP_MA='BLOC',
                          NOM_CMP=('SIXX', 'SIYY', 'SIZZ'),
                          VALE=(-P0, -P0, -P0,), ), );

# ***********************************************************************
#
#    ESSAI TRIAXIAL DRAINE
#
# ***********************************************************************

U1 = STAT_NON_LINE(MODELE=MODELE,
                   CHAM_MATER=CHMAT,
                   EXCIT=(_F(CHARGE=SIGLAT, FONC_MULT=COEF3, ),
                          _F(CHARGE=DEPHAUT, FONC_MULT=COEF2, ),
                          _F(CHARGE=DEPL_1, FONC_MULT=COEF3, ),),
                   ETAT_INIT=_F(SIGM=SIG0, ),
                   COMPORTEMENT=_F(RELATION='RANKINE', ),
                   NEWTON=_F(MATRICE='TANGENTE',
                             PREDICTION='ELASTIQUE',
                             REAC_ITER=1, ),
                   CONVERGENCE=_F(RESI_GLOB_RELA=1.E-6,
                                  ITER_GLOB_MAXI=30,
                                  ARRET='NON', ),
                   SOLVEUR=_F(METHODE='MUMPS', NPREC=8, ),
                   INCREMENT=_F(LIST_INST=TEMPS, ), );

U2 = MECA_NON_LINE(MODELE=MODELE,
                   CHAM_MATER=CHMAT,
                   EXCIT=(_F(CHARGE=SIGLAT, FONC_MULT=COEF3, ),
                          _F(CHARGE=DEPHAUT, FONC_MULT=COEF2, ),
                          _F(CHARGE=DEPL_1, FONC_MULT=COEF3, ),),
                   ETAT_INIT=_F(SIGM=SIG0, ),
                   COMPORTEMENT=_F(RELATION='RANKINE', ),
                   NEWTON=_F(MATRICE='TANGENTE',
                             PREDICTION='ELASTIQUE',
                             REAC_ITER=1, ),
                   CONVERGENCE=_F(RESI_GLOB_RELA=1.E-6,
                                  ITER_GLOB_MAXI=30,
                                  ARRET='NON', ),
                   SOLVEUR=_F(METHODE='MUMPS', NPREC=8, ),
                   INCREMENT=_F(LIST_INST=TEMPS, ), );

# TEST ETAT_INIT  AVEC SIGM
# *************************

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant

nbRank = U1.getNumberOfRanks()

DEPL_REF = U1.getFieldOnNodesReal("DEPL", nbRank - 1)
SIGMA_REF = U1.getFieldOnCellsReal("SIEF_ELGA", nbRank - 1)
VARI_REF = U1.getFieldOnCellsReal("VARI_ELGA", nbRank - 1)

DEPL = U2.getFieldOnNodesReal("DEPL", nbRank - 1)
SIGMA = U2.getFieldOnCellsReal("SIEF_ELGA", nbRank - 1)
VARI = U2.getFieldOnCellsReal("VARI_ELGA", nbRank - 1)

DIF_DEPL = DEPL_REF - DEPL
DIF_SIG = SIGMA_REF - SIGMA
DIF_VAR = VARI_REF - VARI

TEST_RESU(CHAM_ELEM=(_F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        PRECISION=1e-5,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_SIG,
                        VALE_CALC=1E-10,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        PRECISION=1e-5,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_SIG,
                        VALE_CALC=1E-10,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        ORDRE_GRANDEUR=10e-3,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_VAR,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        ORDRE_GRANDEUR=10e-3,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_VAR,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     ),
          )

# EVOL_NOLI
# *********

U3 = STAT_NON_LINE(MODELE=MODELE,
                   CHAM_MATER=CHMAT,
                   EXCIT=(_F(CHARGE=SIGLAT, FONC_MULT=COEF3, ),
                          _F(CHARGE=DEPHAUT, FONC_MULT=COEF2, ),
                          _F(CHARGE=DEPL_1, FONC_MULT=COEF3, ),),
                   ETAT_INIT=_F(EVOL_NOLI=U1, INST_ETAT_INIT=14),
                   COMPORTEMENT=_F(RELATION='RANKINE', ),
                   NEWTON=_F(MATRICE='TANGENTE',
                             PREDICTION='ELASTIQUE',
                             REAC_ITER=1, ),
                   CONVERGENCE=_F(RESI_GLOB_RELA=1.E-6,
                                  ITER_GLOB_MAXI=10,
                                  ARRET='NON', ),
                   SOLVEUR=_F(METHODE='MUMPS', NPREC=8, ),
                   INCREMENT=_F(LIST_INST=TEMPS2, ), );

U4 = MECA_NON_LINE(MODELE=MODELE,
                   CHAM_MATER=CHMAT,
                   EXCIT=(_F(CHARGE=SIGLAT, FONC_MULT=COEF3, ),
                          _F(CHARGE=DEPHAUT, FONC_MULT=COEF2, ),
                          _F(CHARGE=DEPL_1, FONC_MULT=COEF3, ),),
                   ETAT_INIT=_F(EVOL_NOLI=U2, INST_ETAT_INIT=14),
                   COMPORTEMENT=_F(RELATION='RANKINE', ),
                   NEWTON=_F(MATRICE='TANGENTE',
                             PREDICTION='ELASTIQUE',
                             REAC_ITER=1, ),
                   CONVERGENCE=_F(RESI_GLOB_RELA=1.E-8,
                                  ITER_GLOB_MAXI=10,
                                  ARRET='NON', ),
                   SOLVEUR=_F(METHODE='MUMPS', NPREC=8, ),
                   INCREMENT=_F(LIST_INST=TEMPS2, ), );

# TEST EVOL_NOLI et INST_ETAT_INIT
# *******************************

# Test des intervalles de temps
test = code_aster.TestCase()

nbRank3 = U3.getNumberOfRanks()
nbRank4 = U4.getNumberOfRanks()

range3 = [U3.getTimeValue(i) for i in range(nbRank3)]
range4 = [U3.getTimeValue(i) for i in range(nbRank4)]

test.assertEqual(nbRank3, nbRank4)
test.assertEqual(range3, range4)

DEPL_REF = U3.getFieldOnNodesReal("DEPL", nbRank3 - 1)
SIGMA_REF = U3.getFieldOnCellsReal("SIEF_ELGA", nbRank3 - 1)
VARI_REF = U3.getFieldOnCellsReal("VARI_ELGA", nbRank3 - 1)

DEPL = U4.getFieldOnNodesReal("DEPL", nbRank3 - 1)
SIGMA = U4.getFieldOnCellsReal("SIEF_ELGA", nbRank3 - 1)
VARI = U4.getFieldOnCellsReal("VARI_ELGA", nbRank3 - 1)

DIF_DEPL = DEPL_REF - DEPL
DIF_SIG = SIGMA_REF - SIGMA
DIF_VAR = VARI_REF - VARI

TEST_RESU(CHAM_ELEM=(_F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        PRECISION=1e-5,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_SIG,
                        VALE_CALC=1E-10,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        PRECISION=1e-5,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_SIG,
                        VALE_CALC=1E-10,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        ORDRE_GRANDEUR=10e-3,
                        TYPE_TEST='MIN',
                        CHAM_GD=DIF_VAR,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     _F(CRITERE='ABSOLU',
                        REFERENCE='ANALYTIQUE',
                        ORDRE_GRANDEUR=10e-3,
                        TYPE_TEST='MAX',
                        CHAM_GD=DIF_VAR,
                        VALE_CALC=0.0,
                        VALE_REFE=0.0,
                        VALE_ABS='OUI', ),
                     ),
          )

TEST_RESU(CHAM_NO=(_F(CRITERE='ABSOLU',
                      REFERENCE='ANALYTIQUE',
                      PRECISION=1e-8,
                      TYPE_TEST='MIN',
                      CHAM_GD=DIF_DEPL,
                      VALE_CALC=1.e-10,
                      VALE_REFE=0.0,
                      VALE_ABS='OUI', ),
                   _F(CRITERE='ABSOLU',
                      REFERENCE='ANALYTIQUE',
                      PRECISION=1e-8,
                      TYPE_TEST='MAX',
                      CHAM_GD=DIF_DEPL,
                      VALE_CALC=1.e-10,
                      VALE_REFE=0.0,
                      VALE_ABS='OUI', ),
                   ))
FIN();
