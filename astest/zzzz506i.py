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
from code_aster.Solvers import NonLinearSolver, TimeStepper
from code_aster import CA
from code_aster.CA import MPI

DEBUT(CODE=_F(NIV_PUB_WEB="INTERNET"), DEBUG=_F(SDVERI="OUI"), INFO=1)

rank = MPI.ASTER_COMM_WORLD.Get_rank()

mesh = LIRE_MAILLAGE(FORMAT="MED", PARTITIONNEUR="PTSCOTCH")
mesh = mesh.refine()

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

FTRACTUB = DEFI_FONCTION(
    NOM_PARA="EPSI",
    VALE=(1000.0, 2000.0, 2000.0, 5000.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))


encast = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=(_F(GROUP_MA="BAS", DX=0, DY=0.0, DZ=0.0),))

depl = AFFE_CHAR_CINE(MODELE=model, MECA_IMPO=(_F(GROUP_MA="HAUT", DZ=1),))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=2))

DEFLIST = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=LIST), ECHEC=_F(ACTION="ARRET"))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

# STAT_NON_LINE DE REFERENCE
myOptions = "-snes_rtol 1.e-14 -ksp_type fgmres  -pc_type lu -ksp_monitor -pc_factor_mat_solver_type mumps -snes_view"
SOLUT = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="GDEF_LOG"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-16),
    METHODE="SNES",
    INCREMENT=_F(LIST_INST=DEFLIST),
    SOLVEUR=_F(METHODE="PETSC", OPTION_PETSC=myOptions),
    INFO=1,
)
solt = SOLUT.getField("DEPL", 1)
solt.printMedFile(f"/tmp/solt_{rank}.resu.med")

myOptions = (
    "-ksp_type fgmres  -pc_type lu -ksp_monitor -pc_factor_mat_solver_type mumps "
    + "-prefix_push lsnes_ "
    + "-snes_monitor -snes_linesearch_type basic -snes_rtol 1.e-8 -snes_atol 1.e-8 -snes_stol 1.e-24 -snes_max_it 20  "
    + "-ksp_type preonly  -pc_type lu -pc_factor_mat_solver_type mumps "
    + "-prefix_pop "
    + "-prefix_push gsnes_  "
    + "-snes_linesearch_type bt -ksp_monitor -ksp_rtol 1.e-8 -ksp_atol 1.e-50  "
    + "-prefix_pop "
)

SOLUN = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE", DEFORMATION="GDEF_LOG"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    METHODE="RASPEN",
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-10),
    SOLVEUR=_F(METHODE="PETSC", OPTION_PETSC=myOptions),
    INCREMENT=_F(LIST_INST=DEFLIST),
    INFO=1,
)
soln = SOLUN.getField("DEPL", 1)

# =========================================================
#          DETERMINATION DE LA REFERENCE
# =========================================================

# ON EXTRAIT LES CHAMPS A TESTER au dernier instant
SIGMA_REF = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_SIEF_R", NOM_CHAM="SIEF_ELGA", RESULTAT=SOLUT, INST=1.0
)

VARI_REF = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_VARI_R", NOM_CHAM="VARI_ELGA", RESULTAT=SOLUT, INST=1.0
)


SIGMA = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_SIEF_R", NOM_CHAM="SIEF_ELGA", RESULTAT=SOLUN, INST=1.0
)

VARI = CREA_CHAMP(
    OPERATION="EXTR", TYPE_CHAM="ELGA_VARI_R", NOM_CHAM="VARI_ELGA", RESULTAT=SOLUN, INST=1.0
)

# =========================================================
#            REALISATION DES TESTS
# =========================================================

DIF_U = soln - solt

DIF_SIG = SIGMA_REF - SIGMA

DIF_VAR = VARI_REF - VARI

TEST_RESU(
    CHAM_NO=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-06,
            ORDRE_GRANDEUR=1.0,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_U,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    ),
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-03,
            ORDRE_GRANDEUR=1.0,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_SIG,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            ORDRE_GRANDEUR=1.0,
            PRECISION=1.0e-03,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    ),
)

FIN()
