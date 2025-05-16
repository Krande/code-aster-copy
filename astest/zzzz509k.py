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

# Copy of hsnv140b test-case with PETIT (and NOT PETIT_REAC) to être plus rapide


from code_aster.Commands import *
from code_aster import CA

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

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


# Annealing
Tdebut = 600.0  # Temperature de debut de restauration
Tfin = 1200.0  # Temperature de restauration complete

# Valeurs du coefficient alpha en fonction de la temperature
annealingCoefFunc = [
    680.0,
    1.154,
    750.0,
    1.335,
    825.0,
    1.308,
    900.0,
    1.165,
    1000.0,
    0.181,
    1075.0,
    0.687,
    1150.0,
    0.408,
]

# Valeurs du coefficient tau_infini en fonction de la temperature
annealingTauFunc = [
    Tdebut,
    1.0,
    680.0,
    0.877,
    750.0,
    0.767,
    825.0,
    0.682,
    900.0,
    0.576,
    1000.0,
    0.071,
    1075.0,
    0.052,
    1150.0,
    0.045,
    Tfin,
    0.0,
]
# ------------------------------------------------
annealingCoef = DEFI_FONCTION(
    NOM_PARA="TEMP", VALE=(annealingCoefFunc), PROL_DROITE="CONSTANT", PROL_GAUCHE="LINEAIRE"
)
#
annealingTau = DEFI_FONCTION(
    NOM_PARA="TEMP", VALE=(annealingTauFunc), PROL_DROITE="CONSTANT", PROL_GAUCHE="CONSTANT"
)
steel = DEFI_MATERIAU(
    ELAS_FO=_F(E=elasticityModulus, NU=PoissonCoef, ALPHA=thermalCoef, TEMP_DEF_ALPHA=20.0),
    ECRO_LINE_FO=_F(D_SIGM_EPSI=hardeningCoef, SY=elasticityYield),
    REST_ECRO=_F(COEF_ECRO=annealingCoef, TAU_INF=annealingTau, TEMP_MINI=Tdebut, TEMP_MAXI=Tfin),
)

timeList = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=800.0, NOMBRE=8))
timeStepper = DEFI_LIST_INST(METHODE="AUTO", DEFI_LIST=_F(LIST_INST=timeList, NB_PAS_MAXI=3))

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
    AFFE=_F(NOM_CHAM="TEMP", LIST_INST=timeList, CHAM_GD=tempField),
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
nonlinearResult = MECA_NON_LINE(
    MODELE=modelMeca,
    CHAM_MATER=materialField,
    EXCIT=_F(CHARGE=clampBC),
    COMPORTEMENT=_F(
        RELATION="VMIS_ISOT_LINE", DEFORMATION="PETIT", TOUT="OUI", POST_INCR="REST_ECRO"
    ),
    INCREMENT=_F(LIST_INST=timeStepper),
    NEWTON=_F(REAC_INCR=1, MATRICE="TANGENTE", REAC_ITER=5),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=50),
)
test.assertEqual(nonlinearResult.getNumberOfIndexes() - 1, 3)

nonlinearResult = MECA_NON_LINE(
    reuse=nonlinearResult,
    MODELE=modelMeca,
    CHAM_MATER=materialField,
    ETAT_INIT=_F(EVOL_NOLI=nonlinearResult),
    EXCIT=_F(CHARGE=clampBC),
    COMPORTEMENT=_F(
        RELATION="VMIS_ISOT_LINE", DEFORMATION="PETIT", TOUT="OUI", POST_INCR="REST_ECRO"
    ),
    INCREMENT=_F(LIST_INST=timeStepper),
    NEWTON=_F(REAC_INCR=1, MATRICE="TANGENTE", REAC_ITER=5),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=50),
)
test.assertEqual(nonlinearResult.getNumberOfIndexes() - 1, 3 * 2)

nonlinearResult = STAT_NON_LINE(
    reuse=nonlinearResult,
    MODELE=modelMeca,
    CHAM_MATER=materialField,
    ETAT_INIT=_F(EVOL_NOLI=nonlinearResult),
    EXCIT=_F(CHARGE=clampBC),
    COMPORTEMENT=_F(
        RELATION="VMIS_ISOT_LINE", DEFORMATION="PETIT", TOUT="OUI", POST_INCR="REST_ECRO"
    ),
    INCREMENT=_F(LIST_INST=timeStepper),
    NEWTON=_F(REAC_INCR=1, MATRICE="TANGENTE", REAC_ITER=5),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=50),
)
test.assertEqual(nonlinearResult.getNumberOfIndexes() - 1, 8)

# Post traitement
nonlinearResult = CALC_CHAMP(
    reuse=nonlinearResult,
    MODELE=modelMeca,
    CHAM_MATER=materialField,
    CONTRAINTE=("SIEF_NOEU"),
    VARI_INTERNE=("VARI_NOEU"),
    RESULTAT=nonlinearResult,
)

# Mêmes valeurs que la version STAT_NON_LINE
variRefe = 6.514218250410631e-06
sigmRefe = 378658087.77723074
TEST_RESU(
    RESU=(
        _F(
            INST=100.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="VARI_NOEU",
            GROUP_NO="NO1",
            NOM_CMP="V1",
            VALE_CALC=6.514218250410631e-06,
            VALE_REFE=variRefe,
            REFERENCE="SOURCE_EXTERNE",
            PRECISION=0.0001,
        ),
    )
)
TEST_RESU(
    RESU=(
        _F(
            INST=800.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="SIEF_NOEU",
            GROUP_NO="NO1",
            NOM_CMP="SIYY",
            VALE_CALC=378658087.77723074,
            VALE_REFE=sigmRefe,
            REFERENCE="SOURCE_EXTERNE",
            PRECISION=0.0001,
        ),
    )
)


FIN()
