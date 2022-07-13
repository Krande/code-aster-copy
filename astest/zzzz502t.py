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
import code_aster.LinearAlgebra
from code_aster.Commands import *
from code_aster.Utilities import ExecutionParameter, Options

test = code_aster.TestCase()

code_aster.init("--test")
code_aster.LinearAlgebra.petscInitialize()

# modeling = "D_PLAN_INCO_UPG"
modeling = "D_PLAN_INCO_UP"
part = "PTSCOTCH"  #'SANS'

# -----------------------------------------------------------------------------
# ----------------------------------- Modèle  ---------------------------------
# -----------------------------------------------------------------------------

Mesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR=part)

Model = AFFE_MODELE(
    MAILLAGE=Mesh,
    AFFE=_F(MODELISATION=modeling, PHENOMENE="MECANIQUE", TOUT="OUI"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)

CharCin = AFFE_CHAR_CINE(
    MODELE=Model,
    MECA_IMPO=(
        _F(GROUP_NO="N2", DX=0.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="N4", DX=0.0, DY=0.0, DZ=0.0),
    ),
)

CharMeca = AFFE_CHAR_MECA(
    MODELE=Model,
    LIAISON_DDL=_F(GROUP_NO=("N1", "N3"), DDL=("PRES", "PRES"), COEF_MULT=(1.0, 1.0), COEF_IMPO=0),
    DDL_IMPO=(_F(GROUP_NO="N1", PRES=200000.0),),
    INFO=2,
)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.499999))

AffMat = AFFE_MATERIAU(MAILLAGE=Mesh, AFFE=_F(TOUT="OUI", MATER=MATER1))

# -----------------------------------------------------------------------------
# --------------------------------- Assemblage --------------------------------
# -----------------------------------------------------------------------------

ExecutionParameter().disable(Options.UseLegacyMode)

AssemblyObj = ASSEMBLAGE(
    MODELE=Model,
    CHAM_MATER=AffMat,
    CHARGE=CharMeca,
    CHAR_CINE=CharCin,
    NUME_DDL=CO("asterNume"),
    MATR_ASSE=(_F(MATRICE=CO("asterRigi"), OPTION="RIGI_MECA"),),
)

petscMat = AssemblyObj.asterRigi.toPetsc()
print("Norm: ", petscMat.getSizes())
ref = 1823496.3881588143
test.assertAlmostEqual(petscMat.norm(), ref, delta=ref * 1.0e-6)
test.assertSequenceEqual(petscMat.getSizes(), ((24, 48), (24, 48)))

# petscMat.view()

# -----------------------------------------------------------------------------
# --------------------------------- Fin Aster ---------------------------------
# -----------------------------------------------------------------------------

FIN()
