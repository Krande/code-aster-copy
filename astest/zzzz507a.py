# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
from code_aster import CA

CA.init("--test")
test = CA.TestCase()

# Prepare mesh
mesh = CA.Mesh.buildCube(refine=1)
model = DEFI_GROUP(
    reuse=mesh,
    MAILLAGE=mesh,
    CREA_GROUP_MA=(
        _F(
            OPTION="BANDE",
            NOM="VOLU_HAUT",
            POINT=(0, 0, 0.9),
            VECT_NORMALE=(0.0, 0.0, 1.0),
            DIST=0.15,
        ),
        _F(
            OPTION="BANDE",
            NOM="VOLU_BAS",
            POINT=(0, 0, 0.1),
            VECT_NORMALE=(0.0, 0.0, -1.0),
            DIST=0.15,
        ),
    ),
    INFO=2,
)


model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))
material = DEFI_MATERIAU(ELAS=_F(E=1, NU=0.0, RHO=1))

fieldMate = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=material))

clamp = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="BOTTOM", DX=0.0, DY=0.0, DZ=0.0),))

gravity = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=9.81, DIRECTION=(0.0, 0.0, -1.0)))

numeDof = NUME_DDL(MODELE=model, CHARGE=clamp)

loadElem = CALC_VECT_ELEM(OPTION="CHAR_MECA", CHAM_MATER=fieldMate, CHARGE=gravity)

loadAsse = ASSE_VECTEUR(VECT_ELEM=loadElem, NUME_DDL=numeDof)


maskField = CREA_CHAMP(
    TYPE_CHAM="ELEM_NEUT_I",
    OPERATION="AFFE",
    MODELE=model,
    PROL_ZERO="OUI",
    INFO=2,
    AFFE=(
        _F(GROUP_MA="VOLU_BAS", NOM_CMP="X1", VALE_I=1),
        _F(GROUP_MA="VOLU_HAUT", NOM_CMP="X1", VALE_I=0),
    ),
)

loadAsse2 = loadElem.assembleWithMask(numeDof, maskField, 0)
nodalField2 = loadAsse2.toSimpleFieldOnNodes()
nodalField2.updateValuePointers()

test.assertAlmostEqual(nodalField2[2, 2], -0.15328124999999998)
test.assertAlmostEqual(nodalField2[23, 2], 0)

loadAsse3 = loadElem.assembleWithMask(numeDof, maskField, 1)
nodalField3 = loadAsse3.toSimpleFieldOnNodes()
nodalField3.updateValuePointers()

test.assertAlmostEqual(nodalField3[22, 2], -0.6131250000000001)
test.assertAlmostEqual(nodalField3[1, 2], 0)


test.printSummary()


# IMPR_RESU(FORMAT='MED', RESU=_F(MAILLAGE = mesh))

FIN()
