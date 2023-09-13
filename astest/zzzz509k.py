# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

# Copy of hsnv140b test-case

import code_aster
from code_aster.Commands import *

DEBUT(CODE=_F(NIV_PUB_WEB="INTERNET"), DEBUG=_F(SDVERI="OUI"), INFO=1)

test = code_aster.TestCase()

mesh = LIRE_MAILLAGE(FORMAT="ASTER")

modelMeca = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="AXIS")
)

# Material parameters
elasticityModulus = DEFI_FONCTION(
    NOM_PARA="TEMP",
    VALE=(
        20.0,
        195600.0e6,
        100.0,
        191200.0e6,
        200.0,
        185700.0e6,
        300.0,
        179600.0e6,
        400.0,
        172600.0e6,
        500.0,
        164500.0e6,
        600.0,
        155000.0e6,
        700.0,
        144100.0e6,
        800.0,
        131400.0e6,
        900.0,
        116800.0e6,
        1000.0,
        100000.0e6,
        1100.0,
        80000.0e6,
        1200.0,
        57000.0e6,
        1300.0,
        30000.0e6,
        1400.0,
        2000.0e6,
        1500.0,
        1000.0e6,
    ),
    PROL_DROITE="CONSTANT",
    PROL_GAUCHE="LINEAIRE",
)

PoissonCoef = DEFI_CONSTANTE(VALE=0.3)

thermalCoef = DEFI_FONCTION(
    NOM_PARA="TEMP",
    VALE=(
        20.0,
        14.56e-6,
        100.0,
        15.39e-6,
        200.0,
        16.21e-6,
        300.0,
        16.86e-6,
        400.0,
        17.37e-6,
        500.0,
        17.78e-6,
        600.0,
        18.12e-6,
        700.0,
        18.43e-6,
        800.0,
        18.72e-6,
        900.0,
        18.99e-6,
        1000.0,
        19.27e-6,
        1100.0,
        19.53e-6,
        1200.0,
        19.79e-6,
        1300.0,
        20.02e-6,
        1600.0,
        20.02e-6,
    ),
    PROL_DROITE="CONSTANT",
    PROL_GAUCHE="CONSTANT",
)

elasticityYield = DEFI_FONCTION(
    NOM_PARA="TEMP",
    VALE=(
        20.0,
        286.0e6,
        200.0,
        212.0e6,
        400.0,
        180.0e6,
        600.0,
        137.0e6,
        800.0,
        139.0e6,
        1000.0,
        70.0e6,
        1100.0,
        35.0e6,
        1200.0,
        16.0e6,
        1300.0,
        10.0e6,
        1500.0,
        10.0e6,
    ),
    PROL_DROITE="CONSTANT",
    PROL_GAUCHE="CONSTANT",
)

hardeningCoef = DEFI_FONCTION(
    NOM_PARA="TEMP",
    VALE=(
        20.0,
        2.400e9,
        700.0,
        2.400e9,
        800.0,
        2.350e9,
        900.0,
        1.500e9,
        1000.0,
        0.800e9,
        1100.0,
        0.725e9,
        1200.0,
        0.150e9,
        1300.0,
        0.010e9,
    ),
    PROL_DROITE="CONSTANT",
    PROL_GAUCHE="LINEAIRE",
)


# REST_ECRO
Tdebut = 600.0  # Temperature de debut de restauration
Tfin = 1000.0  # Temperature de restauration complete

restEcroPara = DEFI_FONCTION(
    NOM_PARA="TEMP", VALE=(Tdebut, 1.0, Tfin, 0.0), PROL_DROITE="CONSTANT", PROL_GAUCHE="CONSTANT"
)

steel = DEFI_MATERIAU(
    ELAS_FO=_F(E=elasticityModulus, NU=PoissonCoef, ALPHA=thermalCoef, TEMP_DEF_ALPHA=20.0),
    ECRO_LINE_FO=_F(D_SIGM_EPSI=hardeningCoef, SY=elasticityYield),
    REST_ECRO=_F(FONC_MULT=restEcroPara),
)

timeList = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=800.0, NOMBRE=800))
timeStepper = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=timeList))

# Create thermal result
tempFunc = DEFI_FONCTION(
    NOM_PARA="INST",
    NOM_RESU="TEMP",
    VALE=(
        # cycle 1
        0.0,
        20.0,
        100,
        1125.0,
        200,
        21.0,
        # cycle 2
        300,
        932.0,
        400,
        22.0,
        # cycle 3
        500,
        685.0,
        600,
        22.0,
        # cycle 4
        700,
        473.0,
        800,
        21.0,
    ),
    PROL_GAUCHE="CONSTANT",
    PROL_DROITE="CONSTANT",
)

tempField = CREA_CHAMP(
    OPERATION="AFFE",
    TYPE_CHAM="NOEU_TEMP_F",
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", NOM_CMP="TEMP", VALE_F=tempFunc),
)

tempResult = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    NOM_CHAM="TEMP",
    AFFE=_F(LIST_INST=timeList, CHAM_GD=tempField),
)

# Boundary conditions
clampBC = AFFE_CHAR_MECA(
    MODELE=modelMeca, DDL_IMPO=(_F(GROUP_NO=("NO1", "NO2", "NO3", "NO4"), DY=0.0))
)

# Affect material parameters
materialField = AFFE_MATERIAU(
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", MATER=steel),
    AFFE_VARC=_F(TOUT="OUI", EVOL=tempResult, NOM_VARC="TEMP", NOM_CHAM="TEMP", VALE_REF=20.0),
)

# Mechanic non-linear
# nonlinearResultRefe = STAT_NON_LINE(
#     MODELE=modelMeca,
#     CHAM_MATER=materialField,
#     EXCIT=_F(CHARGE=clampBC),
#     COMPORTEMENT=_F(
#         RELATION="VMIS_ISOT_LINE", DEFORMATION="PETIT_REAC", TOUT="OUI", POST_INCR="REST_ECRO"
#     ),
#     INCREMENT=_F(LIST_INST=timeStepper, NUME_INST_FIN=800),
#     NEWTON=_F(REAC_INCR=1, MATRICE="TANGENTE", REAC_ITER=5),
#     CONVERGENCE=_F(RESI_GLOB_RELA=5.0e-05, ITER_GLOB_MAXI=50),
# )

nonlinearResult = MECA_NON_LINE(
    MODELE=modelMeca,
    CHAM_MATER=materialField,
    EXCIT=_F(CHARGE=clampBC),
    COMPORTEMENT=_F(
        RELATION="VMIS_ISOT_LINE", DEFORMATION="PETIT_REAC", TOUT="OUI", POST_INCR="REST_ECRO"
    ),
    INCREMENT=_F(LIST_INST=timeStepper),
    NEWTON=_F(REAC_INCR=1, MATRICE="TANGENTE", REAC_ITER=5),
    CONVERGENCE=_F(RESI_GLOB_RELA=5.0e-05, ITER_GLOB_MAXI=50),
)

# Post traitement
nonlinearResult = CALC_CHAMP(
    reuse=nonlinearResult,
    MODELE=modelMeca,
    CHAM_MATER=materialField,
    CONTRAINTE=("SIEF_NOEU"),
    VARI_INTERNE=("VARI_NOEU"),
    RESULTAT=nonlinearResult,
)
#
TEST_RESU(
    RESU=(
        _F(
            INST=89.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="VARI_NOEU",
            GROUP_NO="NO1",
            NOM_CMP="V1",
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            REFERENCE="SOURCE_EXTERNE",
            ORDRE_GRANDEUR=1e-6,
            CRITERE="ABSOLU",
        ),
        _F(
            INST=200.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="SIEF_NOEU",
            GROUP_NO="NO1",
            NOM_CMP="SIYY",
            VALE_CALC=313428794.37,
            VALE_REFE=303.0e06,
            REFERENCE="SOURCE_EXTERNE",
            PRECISION=0.1,
        ),
        _F(
            INST=400.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="SIEF_NOEU",
            GROUP_NO="NO1",
            NOM_CMP="SIYY",
            VALE_CALC=312625600.111,
            VALE_REFE=316.0e06,
            REFERENCE="SOURCE_EXTERNE",
            PRECISION=0.1,
        ),
        _F(
            INST=600.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="SIEF_NOEU",
            GROUP_NO="NO1",
            NOM_CMP="SIYY",
            VALE_CALC=311271496.622,
            VALE_REFE=325.0e06,
            REFERENCE="SOURCE_EXTERNE",
            PRECISION=0.1,
        ),
        _F(
            INST=800.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="SIEF_NOEU",
            GROUP_NO="NO1",
            NOM_CMP="SIYY",
            VALE_CALC=336504214.735,
            VALE_REFE=327.0e06,
            REFERENCE="SOURCE_EXTERNE",
            PRECISION=0.1,
        ),
    )
)

FIN()
