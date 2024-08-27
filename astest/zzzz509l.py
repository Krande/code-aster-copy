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

from code_aster.Utilities import petscInitialize

CA.init("--test", "--abort")

test = CA.TestCase()

petscInitialize("-on_error_abort ")

pMesh2 = CA.ParallelMesh.buildRectangle(lx=3, ly=1, nx=3, ny=1)

pMesh2 = pMesh2.refine(4)

model = AFFE_MODELE(
    MAILLAGE=pMesh2, AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI")
)

char_cin = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=(_F(GROUP_MA=("LEFT", "RIGHT"), DX=0.0, DY=0.0),))

acier = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=0.0, RHO=1.0))

chmat = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=acier))

pesa = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=1.0e-1, DIRECTION=(0.0, -1.0, 0.0)))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1))

DEFLIST = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LIST), ECHEC=_F(ACTION="ARRET"))

myOptions = (
    "-ksp_type fgmres  -pc_type lu -ksp_monitor -pc_factor_mat_solver_type mumps "
    + "-prefix_push lsnes_ "
    + "-snes_monitor -snes_linesearch_type basic -snes_rtol 1.e-8 -snes_atol 1.e-50 -snes_stol 1.e-50 -snes_max_it 5 "
    + "-prefix_pop "
    + "-prefix_push gsnes_  "
    + "-snes_linesearch_type basic -ksp_monitor -ksp_rtol 1.e-8 -ksp_atol 1.e-50  "
    + "-prefix_pop "
)

SOLU2 = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=chmat,
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=pesa, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="ELAS", DEFORMATION="GDEF_LOG"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    METHODE="RASPEN",
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-12, RESI_GLOB_MAXI=1e-15, ITER_GLOB_MAXI=10),
    INCREMENT=_F(LIST_INST=DEFLIST),
    SOLVEUR=_F(METHODE="PETSC", OPTION_PETSC=myOptions),
    INFO=1,
)
sol2 = SOLU2.getField("DEPL", 1)


myOptions = "-snes_linesearch_type basic -ksp_type preonly  -pc_type lu -ksp_monitor -pc_factor_mat_solver_type mumps "
SOLU1 = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=chmat,
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=pesa, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="ELAS", DEFORMATION="GDEF_LOG"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    METHODE="SNES",
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-9, ITER_GLOB_MAXI=10),
    INCREMENT=_F(LIST_INST=DEFLIST),
    SOLVEUR=_F(METHODE="PETSC", OPTION_PETSC=myOptions),
    INFO=1,
)
sol1 = SOLU1.getField("DEPL", 1)

DIFF = sol2 - sol1

TEST_RESU(
    CHAM_NO=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-08,
            ORDRE_GRANDEUR=1.0e-3,
            TYPE_TEST="MAX",
            CHAM_GD=DIFF,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)


FIN()
