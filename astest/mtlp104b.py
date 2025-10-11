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

import math

####################################################################################################
#
#### Paramètres
#
####################################################################################################

# parametres du modele (cf note EC2M RAP-EDF-CQN02549CE2_05/12/09)
bainCoefB = 0.064380728515158
bainCoefN = 0.216028726635953
#
b_martensite = 0.22414091
n_martensite = 0.14349531
#
R = 2.0
delta_H = 100000.0
C0 = (math.log(10) * R) / delta_H

ac3 = 802.60
tempTemperingBegin = 700.0
tempTemperingEnd = 550.0
tempHold = 610.0

####################################################################################################
#
#### Calcul
#
####################################################################################################

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))

# Get mesh
meshMeca = LIRE_MAILLAGE(FORMAT="MED")

# Material
INCLUDE(UNITE=12)

materialField = AFFE_MATERIAU(MAILLAGE=meshMeca, AFFE=(_F(TOUT="OUI", MATER=steel),))

# FE model
thermalModel = AFFE_MODELE(
    MAILLAGE=meshMeca, AFFE=_F(TOUT="OUI", PHENOMENE="THERMIQUE", MODELISATION=("3D",))
)

# Création artificielle d'un résultat thermique
temp0 = 20.0
temp1 = 150
temp2 = tempTemperingEnd
temp3 = tempTemperingBegin

fieldTemp0 = CREA_CHAMP(
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    MAILLAGE=meshMeca,
    AFFE=_F(TOUT="OUI", NOM_CMP=("TEMP"), VALE=(temp0)),
)
fieldTemp1 = CREA_CHAMP(
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    MAILLAGE=meshMeca,
    AFFE=_F(TOUT="OUI", NOM_CMP=("TEMP"), VALE=(temp1,)),
)

fieldTemp2 = CREA_CHAMP(
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    MAILLAGE=meshMeca,
    AFFE=_F(TOUT="OUI", NOM_CMP=("TEMP"), VALE=(temp2)),
)

fieldTemp3 = CREA_CHAMP(
    TYPE_CHAM="NOEU_TEMP_R",
    OPERATION="AFFE",
    MAILLAGE=meshMeca,
    AFFE=_F(TOUT="OUI", NOM_CMP=("TEMP"), VALE=(temp3)),
)

result = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    AFFE=(
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp0,
            INST=0.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp1,
            INST=1.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp2,
            INST=2.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp3,
            INST=3.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
    ),
)

result2 = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    AFFE=(
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp0,
            INST=0.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp1,
            INST=1.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp2,
            INST=2.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
        _F(
            NOM_CHAM="TEMP",
            CHAM_GD=fieldTemp3,
            INST=3.0,
            MODELE=thermalModel,
            CHAM_MATER=materialField,
        ),
    ),
)

L_inst0 = result.LIST_VARI_ACCES()["INST"]

####################################################################################################
#
#### Calcul métallurgique avec revenu
#
####################################################################################################
PHASINIT = CREA_CHAMP(
    TYPE_CHAM="CART_VAR2_R",
    OPERATION="AFFE",
    MAILLAGE=meshMeca,
    AFFE=_F(
        TOUT="OUI",
        NOM_CMP=("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
        VALE=(0.7, 0.0, 0.3, 0.0, 0.0, 1.0, 10.0),
    ),
)


# V1 : FERRITE
# V2 : PERLITE
# V3 : BAINITE
# V4 : BAIN_REVENU
# V5 : MARTENSITE
# V6 : MARTENSITE_REVEN
# V7 : AUSTENITE
# V8 : SUMCOLD
# V9 : TAILLE_GRAIN
# V10 : TEMP
# V11 : TEMP_MARTENSITE
# V12 : CYCL_THER

result = CALC_META(
    reuse=result,
    MODELE=thermalModel,
    CHAM_MATER=materialField,
    RESULTAT=result,
    ETAT_INIT=_F(META_INIT_ELNO=PHASINIT),
    COMPORTEMENT=_F(RELATION="ACIER", TOUT="OUI", LOI_META="WAECKEL"),
    REVENU=_F(RELATION="ACIER_REVENU", TOUT="OUI", LOI_META="JMA"),
    OPTION=("META_ELNO", "DURT_ELNO"),
)

####################################################################################################
#
#### Tests avec revenu
#
####################################################################################################

metaWithTempering = result.getField("META_ELNO", para="INST", value=L_inst0[-1])

ferriteWith = 0.7
bainiteBaseWith = 0.2700333803022067
martensiteBaseWith = 0.0
bainiteTemperWith = 0.029966619697793316
martensiteTemperWith = 0.0
austeniteWith = 0.0

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V1",
            POINT=1,
            VALE_CALC=ferriteWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=1,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=2,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=3,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=4,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=5,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=6,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=7,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V3",
            POINT=8,
            VALE_CALC=bainiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V5",
            POINT=1,
            VALE_CALC=martensiteBaseWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V7",
            POINT=1,
            VALE_CALC=austeniteWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V4",
            POINT=1,
            VALE_CALC=bainiteTemperWith,
        ),
        _F(
            CHAM_GD=metaWithTempering,
            GROUP_MA="CellToTest",
            NOM_CMP="V6",
            POINT=1,
            VALE_CALC=martensiteTemperWith,
        ),
    )
)

hardness = 230.0

hardnessField = result.getField("DURT_ELNO", para="INST", value=L_inst0[-1])
TEST_RESU(
    CHAM_ELEM=(
        _F(CHAM_GD=hardnessField, GROUP_MA="CellToTest", NOM_CMP="HV", POINT=1, VALE_CALC=hardness),
    )
)

FIN()
