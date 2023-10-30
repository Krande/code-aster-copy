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

import code_aster
from code_aster import LinearSolver, NonLinearResult, PhysicalProblem
from code_aster.Commands import *
from code_aster.Solvers import NonLinearSolver, ProblemSolver
from code_aster.Solvers import SolverOptions as SOP
from code_aster.Solvers import TimeStepper, ProblemType
from code_aster.Solvers.StepSolvers import StaticStepSolver

DEBUT(
    CODE=_F(NIV_PUB_WEB="INTERNET"), ERREUR=_F(ALARME="EXCEPTION"), DEBUG=_F(SDVERI="OUI"), INFO=1
)

test = code_aster.TestCase()

mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

FTRACTUB = DEFI_FONCTION(
    NOM_PARA="EPSI",
    VALE=(1000.0, 2000.0, 2000.0, 5000.0),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

# Very high elasticity limit to simulate elasticity
acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200000.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))


encast = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="BAS", DX=0, DY=0.0, DZ=0.0),))

depl = AFFE_CHAR_MECA(MODELE=model, DDL_IMPO=(_F(GROUP_MA="HAUT", DZ=1.0),))

LIST = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=2))

RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1000.0, 1000.0))

# STAT_NON_LINE DE REFERENCE
SOLUT = STAT_NON_LINE(
    MODELE=model,
    CHAM_MATER=mater,
    EXCIT=(_F(CHARGE=encast, FONC_MULT=RAMPE), _F(CHARGE=depl, FONC_MULT=RAMPE)),
    COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    INCREMENT=_F(LIST_INST=LIST),
    INFO=1,
)


class CustomStepSolver(StaticStepSolver):
    """Example of custom object: just add a print.

    Whatever type of inherited from BaseFeature object can be used (not necessarly
    StepSolver inherited).
    It must provide the StepSolver feature and must have the expected interface.
    """

    def solve(self):
        """Solve the step."""
        try:
            super().solve()
        except code_aster.ConvergenceError as exc:
            print(f"+++ CustomStepSolver raises ConvergenceError: {exc}")
            raise
        else:
            print("+++ CustomStepSolver ends successfully, time:", self.phys_state.time_curr)


def post_hook(nl_solver):
    """Example of hook function.

    Arguments:
        phys_state (PhysicalState): Current physical state.
    """
    print(f"calling hook at time = {nl_solver.phys_state.time_curr}...", flush=True)


class PostHook:
    """User object to be used as a PostStepHook."""

    provide = SOP.PostStepHook

    def __call__(self, nl_solver):
        """Example of hook."""
        nl_solver.phys_state.debugPrint()


snl = ProblemSolver(NonLinearSolver(), NonLinearResult(), pb_type=ProblemType.Static)
snl.use(PhysicalProblem(model, mater))
snl.use(LinearSolver.factory(METHODE="MUMPS"))
snl.phys_pb.addLoadFromDict({"CHARGE": encast, "FONC_MULT": RAMPE})
snl.phys_pb.addLoadFromDict({"CHARGE": depl, "FONC_MULT": RAMPE})
snl.setKeywords(
    CONVERGENCE={"RESI_GLOB_MAXI": 1.0e-6, "ITER_GLOB_MAXI": 20},
    NEWTON={"PREDICTION": "ELASTIQUE"},
    COMPORTEMENT={"RELATION": "VMIS_ISOT_LINE"},
    INFO=1,
)
snl.use(CustomStepSolver())
snl.use(TimeStepper([0.5, 1.0]))
snl.use(post_hook, provide=SOP.PostStepHook)
snl.use(PostHook())
snl.run()

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


SIGMA = snl.phys_state.stress
VARI = snl.phys_state.internVar


# =========================================================
#            REALISATION DES TESTS
# =========================================================


DIF_SIG = SIGMA_REF - SIGMA

DIF_VAR = VARI_REF - VARI

TEST_RESU(
    CHAM_ELEM=(
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-08,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1.5916157281026244e-12,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            PRECISION=1.0e-08,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_SIG,
            VALE_CALC=1.1368683772161603e-12,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            ORDRE_GRANDEUR=1.0e-08,
            TYPE_TEST="MIN",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
        _F(
            CRITERE="ABSOLU",
            REFERENCE="AUTRE_ASTER",
            ORDRE_GRANDEUR=1.0e-8,
            TYPE_TEST="MAX",
            CHAM_GD=DIF_VAR,
            VALE_CALC=0.0,
            VALE_REFE=0.0,
            VALE_ABS="OUI",
        ),
    )
)

FIN()
