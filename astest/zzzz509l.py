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
from code_aster import CA

CA.init("--test", "--abort")
test = CA.TestCase()

mult = 1
mesh = CA.ParallelMesh.buildRectangle(lx=3, ly=1, nx=3 * mult, ny=1 * mult)

mesh = mesh.refine(4)

model = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI")
)

char_cin = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=(_F(GROUP_MA=("LEFT", "RIGHT"), DX=0.0, DY=0.0),))

acier = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=0.0, RHO=1.0))

chmat = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))

pesa = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=1.0e-2, DIRECTION=(0.0, -1.0, 0.0)))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1))

DEFLIST = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LIST), ECHEC=_F(ACTION="ARRET"))


# ------------------------------------------------------------------------------------------------------------------------
# Newton solve

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


# ------------------------------------------------------------------------------------------------------------------------
# RASPEN solve (level 1 only)

myOptions = (
    "-ksp_type fgmres  -pc_type lu -ksp_monitor -pc_factor_mat_solver_type mumps "
    # local snes
    + "-prefix_push lsnes_ "
    + "-snes_monitor -snes_linesearch_type basic -snes_rtol 1.e-8 -snes_atol 1.e-50 -snes_stol 1.e-50 -snes_max_it 5 "
    + "-prefix_pop "
    # global snes
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


# # ------------------------------------------------------------------------------------------------------------------------
# # RASPEN solve with coarse problem

InternalRASPENOpts = "-raspen_with_coarse_pb -raspen_coarse_side Left -raspen_coarse_type SubdomainSVD -raspen_nb_sd_singular_vec 5 "

myOptions = (
    InternalRASPENOpts
    + "-ksp_type preonly -pc_type lu  -snes_monitor -snes_divergence_tolerance -1 -pc_factor_mat_solver_type mumps -snes_linesearch_type basic -snes_max_it 1 "
    + "-snes_atol 1e-50 -snes_rtol 1e-8 -snes_stol 1.e-50 -snes_max_it 5 "
    # local snes
    + "-prefix_push lsnes_ "
    + "-snes_linesearch_type basic -snes_rtol 1.e-7 -snes_atol 1.e-50 -snes_stol 1.e-50 -snes_converged_maxits -snes_monitor -snes_max_it 10 -snes_divergence_tolerance -1 "
    + "-ksp_type preonly  -pc_type lu -pc_factor_mat_solver_type mumps "
    + "-prefix_pop "
    # global snes
    + "-prefix_push gsnes_  "
    + "-snes_linesearch_type basic -snes_rtol 1e-7 -snes_atol 1e-50 -ksp_monitor -ksp_gmres_restart 10000 -ksp_max_it 5000 -ksp_converged_maxits -ksp_rtol 1.e-6 -ksp_atol 1.e-16 -snes_max_it 50 "
    + "-prefix_pop "
)


SOLU3 = MECA_NON_LINE(
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
sol3 = SOLU3.getField("DEPL", 1)


# ------------------------------------------------------------------------------------------------------------------------
# RASPEN solve with substructuring

InternalRASPENOpts = "-raspen_with_substructuring -prefix_push sksp_ -ksp_type gmres -pc_type asm -sub_ksp_type hpddm -sub_ksp_hpddm_type gcrodr -ksp_rtol 1.e-8 -ksp_monitor -sub_ksp_monitor -prefix_pop "

myOptions = (
    InternalRASPENOpts
    + "-ksp_type preonly -pc_type lu  -pc_factor_mat_solver_type mumps -snes_linesearch_type basic "
    # local snes
    + "-prefix_push lsnes_ "
    + "-snes_linesearch_type basic -snes_rtol 1.e-7 -snes_atol 1.e-50 -snes_stol 1.e-50 -snes_monitor -snes_max_it 10 -snes_divergence_tolerance -1 "
    + "-ksp_type preonly  -pc_type lu -pc_factor_mat_solver_type mumps "
    + "-prefix_pop "
    # global snes
    + "-prefix_push gsnes_  "
    + "-snes_linesearch_type basic -snes_rtol 1e-7 -snes_atol 1e-50 -ksp_monitor -ksp_type gmres -ksp_gmres_restart 10000 -ksp_max_it 5000 -ksp_converged_maxits -ksp_rtol 1.e-7 -ksp_atol 1.e-16 -snes_max_it 50 "
    + "-prefix_pop "
)

SOLU4 = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=chmat,
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=pesa, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="ELAS", DEFORMATION="GDEF_LOG"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    METHODE="RASPEN",
    # CONVERGENCE=_F(RESI_GLOB_RELA=1e-10, RESI_GLOB_MAXI=1e-15, ITER_GLOB_MAXI=10),
    INCREMENT=_F(LIST_INST=DEFLIST),
    SOLVEUR=_F(METHODE="PETSC", OPTION_PETSC=myOptions),
    INFO=1,
)
sol4 = SOLU4.getField("DEPL", 1)


# ------------------------------------------------------------------------------------------------------------------------
# Validation

diff = sol2 - sol1
test.assertAlmostEqual(diff.norm(), 0.0)

diff = sol3 - sol1
test.assertAlmostEqual(diff.norm(), 0.0)

diff = sol4 - sol1
test.assertAlmostEqual(diff.norm(), 0.0)


FIN()
