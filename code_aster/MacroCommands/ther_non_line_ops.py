# coding: utf-8

# Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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
from ..Commands import THER_NON_LINE2
from ..Objects import (
    ThermalDirichletBC,
    ThermalLoadFunction,
    ThermalLoadReal,
    ThermalResult,
    ParallelThermalLoadFunction,
    ParallelThermalLoadReal,
    PhysicalProblem,
)
from ..Solvers import NonLinearSolver, ProblemSolver, TimeStepper
from ..Utilities import print_stats, force_list


def use_fortran(keywords):
    excluded_keys = ("EVOL_THER_SECH", "OBSERVATION")

    for key in excluded_keys:
        if key in keywords:
            return True

    if keywords["TYPE_CALCUL"] == "TRAN":
        return True

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

    return False


def ther_non_line_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.
    """

    args = _F(args)

    if use_fortran(args):
        return THER_NON_LINE2(**args)

    # python Version
    if "RESULTAT" in args:
        result = args["RESULTAT"]
    else:
        result = ThermalResult()
    solver = ProblemSolver(NonLinearSolver(), result)

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
        SOLVEUR=args["SOLVEUR"],
    )
    solver.setKeywords(**param)

    # Add stepper
    timeStepper = TimeStepper.from_keywords(**args["INCREMENT"][0])
    solver.use(timeStepper)

    # Run computation
    solver.run()
    print_stats()
    return solver.result
