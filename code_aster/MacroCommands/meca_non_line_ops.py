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
    FrictionType,
    MechanicalDirichletBC,
    MechanicalLoadFunction,
    MechanicalLoadReal,
    NonLinearResult,
    ParallelMechanicalLoadFunction,
    ParallelMechanicalLoadReal,
    PhysicalProblem,
)
from ..Solvers import ContactManager, NonLinearSolver, ProblemSolver, TimeStepper
from ..Utilities import print_stats


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

    if "INCREMENT" in keywords:
        if "NUME_INST_INIT" in keywords["INCREMENT"] or "NUME_INST_FIN" in keywords["INCREMENT"]:
            raise RuntimeError("unsupported value in INCREMENT")
    # FIXME todo: check consistency between INST_INIT and INST_ETAT_INIT

    if "CONVERGENCE" in keywords:
        for key in keywords["CONVERGENCE"]:
            if key in ("RESI_REFE_RELA", "RESI_COMP_RELA"):
                raise RuntimeError("unsupported value in CONVERGENCE: %s" % key)

    if keywords["METHODE"] not in ["NEWTON", "SNES"]:
        raise RuntimeError("unsupported value in METHODE")


def meca_non_line_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.
    """
    UTMESS("A", "QUALITY1_2")

    args = _F(args)

    # Add controls to prohibit unconverted features
    _contact_check(args["CONTACT"])
    _keywords_check(args)
    adapt_for_mgis_behaviour(self, args)

    solver = ProblemSolver(NonLinearSolver(), NonLinearResult())

    phys_pb = PhysicalProblem(args["MODELE"], args["CHAM_MATER"], args["CARA_ELEM"])
    solver.use(phys_pb)

    # Add parameters
    param = dict(
        ARCHIVAGE=args["ARCHIVAGE"],
        COMPORTEMENT=args["COMPORTEMENT"],
        CONTACT=args["CONTACT"],
        CONVERGENCE=args["CONVERGENCE"],
        ETAT_INIT=args["ETAT_INIT"],
        INFO=args["INFO"],
        METHODE=args["METHODE"],
        NEWTON=args["NEWTON"],
        RECH_LINEAIRE=args["RECH_LINEAIRE"],
        SOLVEUR=args["SOLVEUR"],
    )
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

    # Add stepper
    timeStepper = TimeStepper.from_keywords(**args["INCREMENT"][0])
    solver.use(timeStepper)

    # Run computation
    solver.run()
    print_stats()
    return solver.result
