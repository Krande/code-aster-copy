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
    LinearSolver,
    MechanicalDirichletBC,
    MechanicalLoadFunction,
    MechanicalLoadReal,
    NonLinearResult,
    ParallelMechanicalLoadFunction,
    ParallelMechanicalLoadReal,
    PhysicalProblem,
)
from ..Utilities import print_stats
from .NonLinearSolver import (
    ConvergenceManager,
    GeometricSolver,
    IncrementalSolver,
    NonLinearSolver,
    PhysicalState,
    SNESSolver,
    StepSolver,
    StorageManager,
    TimeStepper,
)
from .NonLinearSolver.contact_manager import ContactManager


def _contact_check(CONTACT):
    """Add check to forbid not convered functionnalities"""
    if CONTACT is not None:
        defi = CONTACT["DEFINITION"]

        for zone in defi.getContactZones():
            assert not zone.hasSmoothing
            assert zone.getPairingParameter().getDistanceFunction() is None
            assert zone.getPairingParameter().getElementaryCharacteristics() is None

            if zone.hasFriction:
                assert zone.getFrictionParameter().getType() == FrictionType.Without

        if defi.hasFriction:
            assert CONTACT["ALGO_RESO_FROT"] == "NEWTON"


def _keywords_check(keywords):
    """To forbid unsupported keywords."""

    if "EXCIT" in keywords:
        for load in keywords["EXCIT"]:
            if load["TYPE_CHARGE"] != "FIXE_CSTE":
                raise RuntimeError("TYPE_CHARGE not supported")

    if "INCREMENT" in keywords:
        if "NUME_INST_INIT" in keywords["INCREMENT"] or "NUME_INST_FIN" in keywords["INCREMENT"]:
            raise RuntimeError("unsupported value in INCREMENT")
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

    # add contact check:
    _contact_check(args["CONTACT"])
    _keywords_check(args)

    snl = NonLinearSolver()
    phys_state = PhysicalState()
    snl.use(phys_state)

    phys_pb = PhysicalProblem(args["MODELE"], args["CHAM_MATER"], args["CARA_ELEM"])
    snl.use(phys_pb)

    snl.use(StorageManager(NonLinearResult()))

    snl.setLoggingLevel(args["INFO"])

    # Add parameters
    param = dict(
        CONVERGENCE=args["CONVERGENCE"],
        NEWTON=args["NEWTON"],
        ETAT_INIT=args["ETAT_INIT"],
        INCREMENT=args["INCREMENT"],
        INFO=args["INFO"],
        CONTACT=args["CONTACT"],
        METHODE=args["METHODE"],
    )
    snl.setKeywords(**param)

    # Add behaviour
    adapt_for_mgis_behaviour(self, args)
    snl.setBehaviourProperty(args["COMPORTEMENT"])

    # Add loads
    if args["EXCIT"] is not None:
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
        definition = args["CONTACT"]["DEFINITION"]
        contact_manager = ContactManager(definition, phys_pb)
        fed_defi = definition.getFiniteElementDescriptor()
        phys_pb.getListOfLoads().addContactLoadDescriptor(fed_defi, None)

    # Add stepper
    timeStepper = TimeStepper(args["INCREMENT"]["LIST_INST"].getValues()[1::])
    if "INST_INIT" in args["INCREMENT"]:
        timeStepper.setInitialStep(args["INCREMENT"]["INST_INIT"], args["INCREMENT"]["PRECISION"])

    if "INST_FIN" in args["INCREMENT"]:
        timeStepper.setFinalStep(args["INCREMENT"]["INST_FIN"], args["INCREMENT"]["PRECISION"])

    if args["ETAT_INIT"] is not None:
        if "EVOL_NOLI" in args["ETAT_INIT"]:
            resu = args["ETAT_INIT"].get("EVOL_NOLI")
            assert isinstance(resu, NonLinearResult), resu
            tini = resu.getTimeValue(resu.getNumberOfIndexes() - 1)
            if "INST_ETAT_INIT" in args["ETAT_INIT"]:
                tini = args["ETAT_INIT"].get("INST_ETAT_INIT")
            timeStepper.setInitialStep(tini)

    snl.use(timeStepper)

    # Define convergence object for the incremental solver
    conv_crit = ConvergenceManager()
    conv_crit.use(phys_pb)
    conv_crit.use(phys_state)
    for crit in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI"):
        epsilon = args["CONVERGENCE"].get(crit)
        if epsilon is not None:
            conv_crit.addCriteria(crit, epsilon)
    if args["CONTACT"]:
        if args["CONTACT"].get("ALGO_RESO_GEOM") == "NEWTON":
            conv_crit.addCriteria("RESI_GEOM", args["CONTACT"].get("RESI_GEOM"))

    # Define the linear solver
    linear_solver = LinearSolver.factory("STAT_NON_LINE", args["SOLVEUR"])

    # Define incremental solver
    incr_solver = IncrementalSolver()
    incr_solver.use(phys_pb)
    incr_solver.use(phys_state)
    incr_solver.use(linear_solver)
    incr_solver.use(contact_manager)
    incr_solver.use(conv_crit)

    # Define the step convergence criteria
    if args["METHODE"] == "NEWTON":
        step_conv_solv = GeometricSolver()
    else:
        step_conv_solv = SNESSolver()
    step_conv_solv.use(phys_pb)
    step_conv_solv.use(phys_state)
    step_conv_solv.setParameters(param)
    step_conv_solv.use(contact_manager)
    step_conv_solv.use(conv_crit)
    step_conv_solv.use(incr_solver)

    # Define convergence object for the step solver
    step_crit = ConvergenceManager()
    step_crit.use(phys_pb)
    step_crit.use(phys_state)
    epsilon = 1.0e150
    epsilon = (args["CONTACT"] or {}).get("RESI_GEOM", epsilon)
    step_crit.addCriteria("RESI_GEOM", epsilon)

    # Define the step solver
    step_solver = StepSolver()
    step_solver.use(phys_pb)
    step_solver.use(phys_state)
    step_solver.use(step_conv_solv)
    step_solver.setParameters(param)
    step_solver.use(step_crit)
    snl.use(step_solver)

    # Run computation
    snl.run()
    print_stats()
    return snl.getResult()
