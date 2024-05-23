# coding: utf-8

# Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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
from ..Helpers.syntax_adapters import adapt_increment_init
from ..Messages import UTMESS
from ..Objects import (
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
from ..Solvers.Post import Annealing, ComputeDisplFromHHO
from ..Utilities import print_stats, reset_stats


def _contact_check(CONTACT):
    """Add controls to prohibit unconverted features in contact"""
    if CONTACT:
        assert CONTACT[0]["ALGO_RESO_GEOM"] == "NEWTON"

        defi = CONTACT[0]["DEFINITION"]

        for zone in defi.getContactZones():
            assert not zone.hasSmoothing
            assert zone.getPairingParameter().getDistanceFunction() is None
            assert zone.getPairingParameter().getElementaryCharacteristics() is None

        if defi.hasFriction:
            assert CONTACT[0]["ALGO_RESO_FROT"] == "NEWTON"


def _keywords_check(keywords):
    """Add controls to prohibit unconverted features."""

    if "EXCIT" in keywords:
        for load in keywords["EXCIT"]:
            if load["TYPE_CHARGE"] not in ("FIXE_CSTE", "DIDI"):
                raise RuntimeError("TYPE_CHARGE not supported")

    if "CONVERGENCE" in keywords:
        for key in keywords["CONVERGENCE"]:
            if key in ("RESI_REFE_RELA", "RESI_COMP_RELA"):
                raise RuntimeError("unsupported value in CONVERGENCE: %s" % key)

    if keywords["METHODE"] not in ["NEWTON", "SNES", "RASPEN"]:
        raise RuntimeError("unsupported value in METHODE")


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
        "COMPORTEMENT": args["COMPORTEMENT"],
        "CONTACT": args["CONTACT"],
        "CONVERGENCE": args["CONVERGENCE"][0],
        "ETAT_INIT": args["ETAT_INIT"] and args["ETAT_INIT"][0],
        "INFO": args["INFO"],
        "METHODE": args["METHODE"],
        "NEWTON": args["NEWTON"][0],
        "RECH_LINEAIRE": args["RECH_LINEAIRE"] and args["RECH_LINEAIRE"][0],
        "SOLVEUR": args["SOLVEUR"][0],
        "REUSE": args["reuse"],
        "INCREMENT": args["INCREMENT"][0],
    }

    print('coucou param["NEWTON"]',param["NEWTON"])

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

    # Register hooks
    solver.use(Annealing())
    solver.use(ComputeDisplFromHHO())

    # Run computation
    solver.run()
    print_stats()
    reset_stats()
    return solver.result
