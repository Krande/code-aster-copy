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

# ***********************************************************************
#
# Initial state for elastic behaviours
#
# ***********************************************************************

DEBUT(CODE="OUI")

# ***********************************************************************
#
# Mesh and model
#
# ***********************************************************************

MAILLAGE = LIRE_MAILLAGE(FORMAT="MED", PARTITIONNEUR="PTSCOTCH")

MODELE = AFFE_MODELE(
    MAILLAGE=MAILLAGE, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="AXIS")
)

MAILLAGE = MODI_MAILLAGE(
    reuse=MAILLAGE,
    MAILLAGE=MAILLAGE,
    ORIE_PEAU=_F(GROUP_MA_PEAU=("DROIT", "GAUCHE", "BAS", "HAUT")),
    INFO=1,
)

# ***********************************************************************
#
# Time discretization
#
# ***********************************************************************

tarret = 10.0
tfin = 20.0
npas = 30
temps_max = 30.0

dtemps = temps_max / npas

ltemps = [dtemps * i for i in range(npas + 1)]

TEMPS = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=20, NOMBRE=10),))
TEMPS2 = DEFI_LIST_REEL(DEBUT=20, INTERVALLE=(_F(JUSQU_A=temps_max, NOMBRE=5),))

# ***********************************************************************
#
# Material parameters
#
# ***********************************************************************

# modules mecaniques [kPa]
# ------------------------
# K=516.2E6
# G=238.2E6
# # =>
# YOUNG = 9.*K*G /(3.*K+G)
# POISSON = (3.*K-2.*G) /(6.*K+2.*G)

YOUNG = 1e6
POISSON = 0.25
SIGMA_T = 1.0e3
K = YOUNG / 3.0 / (1.0 - 2.0 * POISSON)
G = YOUNG / 2.0 / (1.0 + POISSON)

SOL = DEFI_MATERIAU(
    ELAS=_F(E=YOUNG, NU=POISSON, ALPHA=0.0),
    ECRO_LINE=_F(D_SIGM_EPSI=1, SY=1e10),
    RANKINE=_F(SIGMA_T=SIGMA_T),
    INFO=1,
)

CHMAT = AFFE_MATERIAU(MAILLAGE=MAILLAGE, AFFE=_F(TOUT="OUI", MATER=SOL))

# ***********************************************************************
#
# Loads and boundary conditions
#
# ***********************************************************************

# pression de preconsolidation [en kPa]
P0 = 5.0e3
EPZZ = 0.03

npas = 300
dtemps = temps_max / npas
linst = [dtemps * i for i in range(npas)]

# SIGLAT = AFFE_CHAR_MECA(MODELE=MODELE, PRES_REP=_F(GROUP_MA=("DROIT",), PRES=P0))

# DEPHAUT = AFFE_CHAR_CINE(MODELE=MODELE, MECA_IMPO=(_F(GROUP_MA=("HAUT",), DY=1.0),))

DEPL_1 = AFFE_CHAR_CINE(
    MODELE=MODELE, MECA_IMPO=(_F(GROUP_MA="BAS", DY=0.0), _F(GROUP_NO="B", DX=0.0))
)


# COEF2 = DEFI_FONCTION(
#     NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 0.0, 20, EPZZ, temps_max, 0)
# )

rampUnit = DEFI_FONCTION(NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 1.0))


# ***********************************************************************
#
#  Initial stress/displacements
#
# ***********************************************************************

SIG0 = CREA_CHAMP(
    INFO=2,
    TYPE_CHAM="ELGA_SIEF_R",
    OPERATION="AFFE",
    MODELE=MODELE,
    PROL_ZERO="OUI",
    AFFE=_F(GROUP_MA="BLOC", NOM_CMP=("SIXX", "SIYY", "SIZZ"), VALE=(-P0, -P0, -P0)),
)


# U0 = CREA_CHAMP(
#     INFO=2,
#     TYPE_CHAM="NOEU_DEPL_R",
#     OPERATION="AFFE",
#     MODELE=MODELE,
#     AFFE=_F(GROUP_MA="BLOC", NOM_CMP=("DX", "DY", "DZ"), VALE=(-0.001, -0.001, -0.001)),
# )
# ***********************************************************************
#
#  Computation with STAT_NON_LINE
#
# ***********************************************************************

U1 = STAT_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(_F(CHARGE=DEPL_1, FONC_MULT=rampUnit),),
    ETAT_INIT=_F(SIGM=SIG0),
    COMPORTEMENT=_F(RELATION="ELAS"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=30),
    INCREMENT=_F(LIST_INST=TEMPS, INST_FIN=2.0),
)

# ***********************************************************************
#
#  Computation with MECA_NON_LINE
#
# ***********************************************************************

U2 = MECA_NON_LINE(
    MODELE=MODELE,
    CHAM_MATER=CHMAT,
    EXCIT=(_F(CHARGE=DEPL_1, FONC_MULT=rampUnit),),
    ETAT_INIT=_F(SIGM=SIG0),
    COMPORTEMENT=_F(RELATION="ELAS"),
    NEWTON=_F(MATRICE="ELASTIQUE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_MAXI=1.0e-6, ITER_GLOB_MAXI=30),
    INCREMENT=_F(LIST_INST=TEMPS, INST_FIN=2.0),
)

DEPL1 = U1.getField("DEPL", 1).getValues()
print(DEPL1)

SIEF1 = U1.getField("SIEF_ELGA", 1).getValues()
print(SIEF1)


DEPL2 = U2.getField("DEPL", 1).getValues()
print(DEPL2)

SIEF2 = U2.getField("SIEF_ELGA", 1).getValues()
print(SIEF2)

SIEF2 = U2.getField("SIEF_ELGA", 0).getValues()
print("inst 0")
print(SIEF2)

# # TEST ETAT_INIT  AVEC SIGM
# # *************************

# # ON EXTRAIT LES CHAMPS A TESTER au dernier instant

nbIndexes = U1.getNumberOfIndexes()

DEPL_REF = U1.getField("DEPL", nbIndexes - 1)
SIGMA_REF = U1.getField("SIEF_ELGA", nbIndexes - 1)
VARI_REF = U1.getField("VARI_ELGA", nbIndexes - 1)

DEPL = U2.getField("DEPL", nbIndexes - 1)
SIGMA = U2.getField("SIEF_ELGA", nbIndexes - 1)
VARI = U2.getField("VARI_ELGA", nbIndexes - 1)

DIF_DEPL = DEPL_REF - DEPL
DIF_SIG = SIGMA_REF - SIGMA
DIF_VAR = VARI_REF - VARI

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            PRECISION=1e-5,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1e-10,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-3,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="ANALYTIQUE",
            ORDRE_GRANDEUR=10e-3,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

FIN()
