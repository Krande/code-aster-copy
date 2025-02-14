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
    LinearSolver,
    MechanicalDirichletBC,
    MechanicalLoadFunction,
    MechanicalLoadReal,
    NonLinearResult,
    ParallelContactNew,
    ParallelFrictionNew,
    ParallelMechanicalLoadFunction,
    ParallelMechanicalLoadReal,
    PhysicalProblem,
)
from ..Solvers import (
    BaseOperators,
    ContactManager,
    Context,
    NonLinearOperator,
    PhysicalState,
    TimeStepper,
)
from ..Solvers import ProblemType as PBT
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


def _keyword_clean(obj):
    """Return obj[0] if exists, return obj if not"""

    if hasattr(obj, "__getitem__"):
        if type(obj) is dict:
            return obj
        else:
            return obj[0]

    else:
        return obj


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

    # FIXME: move keywords cleanup in factory - Keep some keywords
    kwds = {
        "ARCHIVAGE": _keyword_clean(args["ARCHIVAGE"]),
        "COMPORTEMENT": args["COMPORTEMENT"],
        "CONTACT": args["CONTACT"],
        "CONVERGENCE": _keyword_clean(args["CONVERGENCE"]),
        "ETAT_INIT": _keyword_clean(args["ETAT_INIT"]),
        "INFO": args["INFO"],
        "METHODE": args["METHODE"],
        "NEWTON": _keyword_clean(args["NEWTON"]),
        "RECH_LINEAIRE": _keyword_clean(args["RECH_LINEAIRE"]),
        "SOLVEUR": _keyword_clean(args["SOLVEUR"]),
        "REUSE": args["reuse"],
        "INCREMENT": _keyword_clean(args["INCREMENT"]),
        "SCHEMA_TEMPS": args.get("SCHEMA_TEMPS"),
    }
    if kwds["SOLVEUR"]["METHODE"] == "PETSC":
        if kwds["SOLVEUR"]["PRE_COND"] == "LDLT_SP":
            kwds["SOLVEUR"]["REAC_PRECOND"] = 0

    phys_pb = PhysicalProblem(args["MODELE"], args["CHAM_MATER"], args["CARA_ELEM"])
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

    operator = NonLinearOperator.factory(phys_pb, result=args.get("reuse"), **kwds)
    operator.run()

    print_stats()
    reset_stats()
    return operator.result
