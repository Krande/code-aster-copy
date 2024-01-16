# coding: utf-8

# Copyright (C) 1991 - 2024  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.


from ..Cata.Syntax import _F
from ..Objects import (
    HHO,
    ParallelThermalLoadFunction,
    ParallelThermalLoadReal,
    PhysicalProblem,
    PostProcessing,
    ThermalDirichletBC,
    ThermalLoadFunction,
    ThermalLoadReal,
    ThermalResult,
)
from ..Solvers import NonLinearSolver, ProblemSolver
from ..Solvers import ProblemType as PBT
from ..Solvers import SolverOptions as SOP
from ..Solvers import TimeStepper
from ..Utilities import force_list, print_stats, reset_stats
from .Utils.ther_non_line_fort_op import THER_NON_LINE_FORT
from ..Helpers.syntax_adapters import adapt_increment_init


def use_fortran(keywords):
    excluded_keys = ("EVOL_THER_SECH", "OBSERVATION")

    # return False

    for key in excluded_keys:
        if key in keywords:
            return True

    # if keywords["TYPE_CALCUL"] == "TRAN":
    #     return True

    if keywords["METHODE"] in ("MODELE_REDUIT", "NEWTON_KRYLOV"):
        return True

    for comp in force_list(keywords["COMPORTEMENT"]):
        if comp["RELATION"] != "THER_NL":
            return True

    affi = keywords["AFFICHAGE"]
    if (
        "UNITE" in affi
        or "PAS" in affi
        or affi["INFO_RESIDU"] == "OUI"
        or affi["INFO_TEMP"] == "OUI"
    ):
        return True

    for load in keywords["EXCIT"]:
        if isinstance(
            load["CHARGE"],
            (
                ThermalLoadFunction,
                ThermalLoadReal,
                # ParallelThermalLoadFunction,
                # ParallelThermalLoadReal,
            ),
        ):
            if load["CHARGE"].hasLoadResult():
                return True

    return False


def ther_non_line_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.
    """

    args = _F(args)
    reset_stats()

    if use_fortran(args):
        return THER_NON_LINE_FORT(**args)

    adapt_increment_init(args, "EVOL_THER")

    verbosity = args.get("INFO") or 1

    # python Version
    if "RESULTAT" in args:
        result = args["RESULTAT"]
    else:
        result = ThermalResult()
    solver = ProblemSolver(NonLinearSolver(), result, PBT.Thermal)

    phys_pb = PhysicalProblem(args["MODELE"], args["CHAM_MATER"], args["CARA_ELEM"])
    # Add loads
    if args["EXCIT"]:
        for load in args["EXCIT"]:
            if isinstance(
                load["CHARGE"],
                (
                    ThermalLoadFunction,
                    ThermalLoadReal,
                    ParallelThermalLoadFunction,
                    ParallelThermalLoadReal,
                    ThermalDirichletBC,
                ),
            ):
                if "FONC_MULT" in load:
                    if isinstance(load["CHARGE"], ThermalDirichletBC):
                        phys_pb.addDirichletBC(load["CHARGE"], load["FONC_MULT"])
                    else:
                        phys_pb.addLoad(load["CHARGE"], load["FONC_MULT"])
                else:
                    if isinstance(load["CHARGE"], ThermalDirichletBC):
                        phys_pb.addDirichletBC(load["CHARGE"])
                    else:
                        phys_pb.addLoad(load["CHARGE"])
            else:
                raise RuntimeError("Unknown load")

    solver.use(phys_pb)

    # Add parameters
    param = dict(
        ARCHIVAGE=args["ARCHIVAGE"],
        COMPORTEMENT=args["COMPORTEMENT"],
        CONVERGENCE=args["CONVERGENCE"],
        ETAT_INIT=args["ETAT_INIT"],
        INFO=args["INFO"],
        METHODE=args["METHODE"],
        NEWTON=args["NEWTON"],
        RECH_LINEAIRE=args["RECH_LINEAIRE"],
        SOLVEUR=args["SOLVEUR"],
        TYPE_CALCUL=args["TYPE_CALCUL"],
        INCREMENT=args["INCREMENT"],
        REUSE=args["reuse"],
    )

    if "SCHEMA_TEMPS" in args:
        param["SCHEMA_TEMPS"] = args["SCHEMA_TEMPS"]

    solver.setKeywords(**param)

    class PostHookHydr:
        """Hook to compute HYDR_ELGA."""

        provide = SOP.PostStepHook

        def __call__(self, nl_solver):
            """Hook to compute HYDR_ELGA"""

            compor = phys_pb.getBehaviourProperty()

            if compor.hasBehaviour("THER_HYDR"):
                post = PostProcessing(nl_solver.phys_pb)

                phys_state = nl_solver.phys_state

                if phys_state._size == 1:
                    hydr_curr = phys_state.createFieldOnCells(nl_solver.phys_pb, "ELGA", "HYDR_R")
                else:
                    hydr_curr = post.computeHydration(
                        phys_state.primal_prev,
                        phys_state.primal_curr,
                        phys_state.time_prev,
                        phys_state.time_curr,
                        phys_state.getState(-1).auxiliary["HYDR_ELGA"],
                    )

                phys_state.set("HYDR_ELGA", hydr_curr)

    class PostHookHHO:
        """Hook to compute HHO_TEMP."""

        provide = SOP.PostStepHook

        def __call__(self, nl_solver):
            """Hook to compute HHO_TEMP"""

            if nl_solver.phys_pb.getModel().existsHHO():
                hho_field = HHO(nl_solver.phys_pb).projectOnLagrangeSpace(
                    nl_solver.phys_state.primal_curr
                )

                nl_solver.phys_state.set("HHO_TEMP", hho_field)

    solver.use(PostHookHydr())
    solver.use(PostHookHHO())

    # Run computation
    solver.run()
    if verbosity:
        print_stats()
    reset_stats()
    return solver.result
