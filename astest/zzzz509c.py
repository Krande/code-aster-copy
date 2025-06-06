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

# Copy of hsnv121b test-case


from code_aster.Commands import *
from code_aster import CA

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

mesh = LIRE_MAILLAGE(FORMAT="ASTER")

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="AXIS"))

elasticityModulus = DEFI_FONCTION(
    NOM_PARA="TEMP", PROL_DROITE="CONSTANT", VALE=(20.0, 250000.0, 120.0, 200000.0)
)

PoissonCoef = DEFI_CONSTANTE(VALE=0.3)

thermalCoef = DEFI_CONSTANTE(VALE=0.0001)

elasticityYield = DEFI_CONSTANTE(VALE=1000.0)

hardeningCoef = DEFI_FONCTION(
    NOM_PARA="TEMP", PROL_DROITE="CONSTANT", VALE=(20.0, 2500.0, 120.0, 2000.0)
)

steel = DEFI_MATERIAU(
    ELAS_FO=_F(E=elasticityModulus, NU=PoissonCoef, ALPHA=thermalCoef, TEMP_DEF_ALPHA=20.0),
    ECRO_LINE_FO=_F(D_SIGM_EPSI=hardeningCoef, SY=elasticityYield),
)

timeStep = DEFI_LIST_REEL(
    DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=1.0, NOMBRE=1), _F(JUSQU_A=2.0, NOMBRE=20))
)

tractionBCFunc = DEFI_FONCTION(
    NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 0.0, 1.0, 0.0, 2.0, 293.3)
)

tempFunc = DEFI_FONCTION(
    NOM_PARA="INST",
    NOM_RESU="TEMP",
    PROL_DROITE="CONSTANT",
    VALE=(0.0, 20.0, 1.0, 120.0, 2.0, 120.0),
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
    AFFE=_F(NOM_CHAM="TEMP", LIST_INST=timeStep, CHAM_GD=tempField),
)

clampBC = AFFE_CHAR_MECA(
    MODELE=model, DDL_IMPO=(_F(GROUP_NO="NO1", DY=0.0), _F(GROUP_NO="NO2", DY=0.0))
)

tractionBC = AFFE_CHAR_MECA(MODELE=model, FACE_IMPO=_F(GROUP_MA="MA2", DY=1.0))

materialField = AFFE_MATERIAU(
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", MATER=steel),
    AFFE_VARC=_F(TOUT="OUI", EVOL=tempResult, NOM_VARC="TEMP", NOM_CHAM="TEMP", VALE_REF=20.0),
)

nonlinearResult = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=materialField,
    EXCIT=(_F(CHARGE=clampBC), _F(CHARGE=tractionBC, FONC_MULT=tractionBCFunc)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="SIMO_MIEHE"),
    INCREMENT=_F(LIST_INST=timeStep),
    NEWTON=_F(PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_MAXI=1.0e-1, ITER_GLOB_MAXI=50),
    # CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=50),
    # RECH_LINEAIRE=_F(ITER_LINE_MAXI=3),
)

# nonlinearResult = STAT_NON_LINE(
#     MODELE=model,
#     CHAM_MATER=materialField,
#     EXCIT=(_F(CHARGE=clampBC), _F(CHARGE=tractionBC, FONC_MULT=tractionBCFunc)),
#     COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="SIMO_MIEHE"),
#     INCREMENT=_F(LIST_INST=timeStep),
#     NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
#     CONVERGENCE=_F(RESI_GLOB_MAXI=1.0e-1, ITER_GLOB_MAXI=50),
#     # RECH_LINEAIRE=_F(ITER_LINE_MAXI=3),
# )

TEST_RESU(
    RESU=(
        _F(
            INST=2.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="DEPL",
            GROUP_NO="NO3",
            NOM_CMP="DY",
            VALE_CALC=293.3,
            VALE_REFE=303.0,
            REFERENCE="NON_DEFINI",
            PRECISION=3.21e-2,
        ),
        _F(
            INST=2.0,
            RESULTAT=nonlinearResult,
            NOM_CHAM="DEPL",
            GROUP_NO="NO3",
            NOM_CMP="DX",
            VALE_CALC=-106.49562908414596,
            VALE_REFE=-110.0,
            REFERENCE="NON_DEFINI",
            PRECISION=3.19e-2,
        ),
        _F(
            INST=2.0,
            POINT=1,
            RESULTAT=nonlinearResult,
            NOM_CHAM="SIEF_ELGA",
            NOM_CMP="SIYY",
            VALE_CALC=1462.0287566770467,
            VALE_REFE=1453.0,
            REFERENCE="NON_DEFINI",
            PRECISION=1.0e-2,
            GROUP_MA="MA1",
        ),
        _F(
            INST=2.0,
            POINT=1,
            RESULTAT=nonlinearResult,
            NOM_CHAM="VARI_ELGA",
            NOM_CMP="V1",
            VALE_CALC=0.25222906115691024,
            VALE_REFE=0.2475,
            REFERENCE="NON_DEFINI",
            PRECISION=1.92e-2,
            GROUP_MA="MA1",
        ),
    )
)

# Ajout d'un test avec un chargement de force extérieures nulles pour tester le calcul RESI_GLOB_RELA
# avec présence de variables de commande TEMP
nonlinearResult2 = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=materialField,
    EXCIT=(_F(CHARGE=clampBC)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="SIMO_MIEHE"),
    INCREMENT=_F(LIST_INST=timeStep, INST_FIN=1.0),
    NEWTON=_F(PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    # CONVERGENCE=_F(RESI_GLOB_MAXI=1.0e-1, ITER_GLOB_MAXI=50),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=50),
    # RECH_LINEAIRE=_F(ITER_LINE_MAXI=3),
)

# nonlinearResult2 = STAT_NON_LINE(
#     MODELE=model,
#     CHAM_MATER=materialField,
#     EXCIT=(_F(CHARGE=clampBC)),
#     COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="SIMO_MIEHE"),
#     INCREMENT=_F(LIST_INST=timeStep, INST_FIN =1.0),
#     NEWTON=_F(PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
#     #CONVERGENCE=_F(RESI_GLOB_MAXI=1.0e-1, ITER_GLOB_MAXI=50),
#     CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=50),
#     # RECH_LINEAIRE=_F(ITER_LINE_MAXI=3),
# )

TEST_RESU(
    RESU=(
        _F(
            INST=1.0,
            RESULTAT=nonlinearResult2,
            NOM_CHAM="DEPL",
            GROUP_NO="NO3",
            NOM_CMP="DY",
            VALE_CALC=9.762840686151652,
            VALE_REFE=9.762840686151653,
            REFERENCE="NON_DEFINI",
            PRECISION=3.21e-2,
        ),
    )
)
FIN()
#
