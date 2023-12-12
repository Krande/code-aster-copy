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
from ..Helpers import adapt_for_mgis_behaviour
from ..Messages import UTMESS
from ..Objects import (
    HHO,
    FrictionType,
    MechanicalDirichletBC,
    MechanicalLoadFunction,
    MechanicalLoadReal,
    NonLinearResult,
    ParallelMechanicalLoadFunction,
    ParallelMechanicalLoadReal,
    PhysicalProblem,
)
from ..Solvers import ContactManager, NonLinearSolver, ProblemSolver
from ..Solvers import ProblemType as PBT
from ..Solvers import SolverOptions as SOP
from ..Utilities import print_stats, reset_stats
from ..Helpers.syntax_adapters import adapt_increment_init


def _contact_check(CONTACT):
    """Add controls to prohibit unconverted features in contact"""
    if CONTACT:
        assert CONTACT[0]["ALGO_RESO_GEOM"] == "NEWTON"

        defi = CONTACT[0]["DEFINITION"]

        for zone in defi.getContactZones():
            assert not zone.hasSmoothing
            assert zone.getPairingParameter().getDistanceFunction() is None
            assert zone.getPairingParameter().getElementaryCharacteristics() is None

            if zone.hasFriction:
                assert zone.getFrictionParameter().getType() == FrictionType.Without

        if defi.hasFriction:
            assert CONTACT[0]["ALGO_RESO_FROT"] == "NEWTON"


def _keywords_check(keywords):
    """Add controls to prohibit unconverted features."""

    if "EXCIT" in keywords:
        for load in keywords["EXCIT"]:
            if load["TYPE_CHARGE"] != "FIXE_CSTE":
                raise RuntimeError("TYPE_CHARGE not supported")

    if "CONVERGENCE" in keywords:
        for key in keywords["CONVERGENCE"]:
            if key in ("RESI_REFE_RELA", "RESI_COMP_RELA"):
                raise RuntimeError("unsupported value in CONVERGENCE: %s" % key)

    if keywords["METHODE"] not in ["NEWTON", "SNES"]:
        raise RuntimeError("unsupported value in METHODE")

    if "COMPORTEMENT" in keywords:
        if "RELATION" in keywords["COMPORTEMENT"]:
            if keywords["COMPORTEMENT"]["RELATION"] != "ELAS":
                raise RuntimeError("unsupported value in RELATION")
        if "DEFORMATION" in keywords["COMPORTEMENT"]:
            if keywords["COMPORTEMENT"]["DEFORMATION"] != "PETIT":
                raise RuntimeError("unsupported value in DEFORMATION")


def meca_non_line_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.
    """
    UTMESS("A", "QUALITY1_2")
    reset_stats()

    args = _F(args)
    adapt_increment_init(args, "EVOL_NOLI")

    # Add controls to prohibit unconverted features
    _contact_check(args["CONTACT"])
    _keywords_check(args)
    adapt_for_mgis_behaviour(self, args)

    # Add parameters
    param = {
        "ARCHIVAGE": args["ARCHIVAGE"],
        "COMPORTEMENT": args["COMPORTEMENT"],
        "CONTACT": args["CONTACT"],
        "CONVERGENCE": args["CONVERGENCE"],
        "ETAT_INIT": args["ETAT_INIT"],
        "INFO": args["INFO"],
        "METHODE": args["METHODE"],
        "NEWTON": args["NEWTON"],
        "RECH_LINEAIRE": args["RECH_LINEAIRE"],
        "SOLVEUR": args["SOLVEUR"],
        "REUSE": args["reuse"],
        "INCREMENT": args["INCREMENT"],
    }

    if param["SOLVEUR"]["METHODE"] == "PETSC":
        if param["SOLVEUR"]["PRE_COND"] == "LDLT_SP":
            param["SOLVEUR"]["REAC_PRECOND"] = 0

    if "SCHEMA_TEMPS" in args:
        problem_type = PBT.MecaDyna
        param["SCHEMA_TEMPS"] = args["SCHEMA_TEMPS"]
    else:
        problem_type = PBT.MecaStat

    result = args.get("reuse")
    if not result:
        result = NonLinearResult()

    # Create the problem solver
    solver = ProblemSolver(NonLinearSolver(), result, problem_type)

    # Create the physical problem (and use it in problem solver)
    phys_pb = PhysicalProblem(args["MODELE"], args["CHAM_MATER"], args["CARA_ELEM"])
    solver.use(phys_pb)

    solver.setKeywords(**param)

    # Add loads
    if args["EXCIT"]:
        for load in args["EXCIT"]:
            if isinstance(
                load["CHARGE"],
                (
                    MechanicalLoadFunction,
                    MechanicalLoadReal,
                    ParallelMechanicalLoadFunction,
                    ParallelMechanicalLoadReal,
                    MechanicalDirichletBC,
                ),
            ):
                phys_pb.addLoadFromDict(load)
            else:
                raise RuntimeError("Unknown load")

    # Add contact
    contact_manager = None
    if args["CONTACT"]:
        definition = args["CONTACT"][0]["DEFINITION"]
        contact_manager = ContactManager(definition, phys_pb)
        fed_defi = definition.getFiniteElementDescriptor()
        phys_pb.getListOfLoads().addContactLoadDescriptor(fed_defi, None)

    solver.use(contact_manager)

    # Add Hook
    class PostHookHHO:
        """User object to be used as a PostStepHook."""

        provide = SOP.PostStepHook

        def __call__(self, nl_solver):
            """Hook to compute HHO_DEPL"""

            if nl_solver.phys_pb.getModel().existsHHO():
                hho_field = HHO(nl_solver.phys_pb).projectOnLagrangeSpace(
                    nl_solver.phys_state.primal_curr
                )
                storage_manager = nl_solver.get_feature(SOP.Storage)
                storage_manager.storeField(
                    nl_solver.step_rank, hho_field, "HHO_DEPL", time=nl_solver.phys_state.time_curr
                )

    solver.use(PostHookHHO())

    # Run computation
    solver.run()
    print_stats()
    reset_stats()
    return solver.result
