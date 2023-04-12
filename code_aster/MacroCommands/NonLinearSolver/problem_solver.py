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

from ...NonLinear import SolverFeature
from ...NonLinear import SolverOptions as SOP
from ...Objects import LinearSolver
from ...Utilities import logger, no_new_attributes
from . import (
    ConvergenceManager,
    GeometricSolver,
    IncrementalSolver,
    PhysicalState,
    SNESSolver,
    StepSolver,
    StorageManager,
)


class ProblemSolver(SolverFeature):
    """Solver for linear and non linear problem.

    Arguments:
        main (*NonLinearFeature*): Main object.
        result (*misc*): The result object.
    """

    required_features = [
        SOP.PhysicalProblem,
        SOP.Storage,
        SOP.Keywords,
        SOP.TimeStepper,
        SOP.LinearSolver,
        SOP.StepSolver,
        SOP.IncrementalSolver,
        SOP.ConvergenceCriteria,
        SOP.ConvergenceManager,
        SOP.ForStep,
        SOP.ForIncr,
    ]
    optional_features = [SOP.Contact]

    _main = _result = None
    _phys_state = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, main, result) -> None:
        super().__init__()
        self._main = main
        self._result = result
        self._phys_state = PhysicalState()

    @property
    def phys_state(self):
        """PhysicalState: current state."""
        return self._phys_state

    @property
    def result(self):
        """*misc* : the result object."""
        return self._result

    def setKeywords(self, **args):
        """Set parameters from user keywords.

        Arguments:
            args (dict) : user keywords.
        """
        # check requirements or set defaults
        self.use(args, SOP.Keywords)

    def initialize(self):
        """Initialize default objects when required."""
        args = self.get_feature(SOP.Keywords, optional=True)
        # FIXME 'use()' or a global logger
        self._main.setLoggingLevel(args.get("INFO", 1))
        # required to build other default objects
        self._main.use(self.phys_pb)
        self._main.use(self._phys_state)
        self._main.use(args, SOP.Keywords)
        self._main.use(self.get_feature(SOP.TimeStepper))
        self._main.use(self._get_storage())
        self._main.use(self._get_step_solver())
        self.check_features()

    def run(self):
        """Solve the problem.

        Returns:
            *misc*: result object.
        """
        self.initialize()
        self._main.run()
        return self._result

    # methods that build and return objects (allowing dependencies)
    def _get_storage(self):
        logger.debug("+++ get StorageManager")
        store = self.get_feature(SOP.Storage, optional=True)
        if not store:
            store = StorageManager(self._result)
        self.use(store)
        return store

    def _get_linear_solver(self):
        logger.debug("+++ get LinearSolver")
        if not self.has_feature(SOP.LinearSolver):
            args = self.get_feature(SOP.Keywords)
            self.use(LinearSolver.factory("STAT_NON_LINE", args["SOLVEUR"]))
        return self.get_feature(SOP.LinearSolver)

    def _get_contact_manager(self):
        logger.debug("+++ get ContactManager")
        return self.get_feature(SOP.Contact, optional=True)

    def _get_incr_conv(self):
        logger.debug("+++ get ConvergenceManager ForIncr")
        converg = self.get_feature(SOP.ConvergenceManager | SOP.ForIncr, optional=True)
        if not converg:
            args = self.get_feature(SOP.Keywords)
            converg = ConvergenceManager()
            for crit in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI"):
                epsilon = args["CONVERGENCE"].get(crit)
                if epsilon is not None:
                    converg.addCriteria(crit, epsilon)
            if args.get("CONTACT"):
                if args["CONTACT"].get("ALGO_RESO_GEOM") == "NEWTON":
                    converg.addCriteria("RESI_GEOM", args["CONTACT"].get("RESI_GEOM"))
        for feat, required in converg.undefined():
            converg.use(self._get(feat | SOP.ForIncr, required))
        self.use(converg, SOP.ForIncr)
        return converg

    def _get_incremental_solver(self):
        logger.debug("+++ get IncrementalSolver")
        incr_solver = self.get_feature(SOP.IncrementalSolver, optional=True)
        if not incr_solver:
            incr_solver = IncrementalSolver()
        for feat, required in incr_solver.undefined():
            incr_solver.use(self._get(feat | SOP.ForIncr, required))
        self.use(incr_solver)
        return incr_solver

    def _get_step_conv(self):
        logger.debug("+++ get ConvergenceManager ForStep")
        step_crit = self.get_feature(SOP.ConvergenceManager | SOP.ForStep, optional=True)
        if not step_crit:
            args = self.get_feature(SOP.Keywords)
            step_crit = ConvergenceManager()
            epsilon = 1.0e150
            epsilon = (args.get("CONTACT") or {}).get("RESI_GEOM", epsilon)
            step_crit.addCriteria("RESI_GEOM", epsilon)
        for feat, required in step_crit.undefined():
            step_crit.use(self._get(feat, required))
        self.use(step_crit, SOP.ForStep)
        return step_crit

    def _get_step_conv_solver(self):
        logger.debug("+++ get ConvergenceCriteria")
        step_conv_solv = self.get_feature(SOP.ConvergenceCriteria, optional=True)
        args = self.get_feature(SOP.Keywords)
        if not step_conv_solv:
            if args.get("METHODE", "NEWTON") == "NEWTON":
                step_conv_solv = GeometricSolver()
            else:
                step_conv_solv = SNESSolver()
        # FIXME replace by 'use(*, Keywords)'
        if not step_conv_solv.param:
            step_conv_solv.setParameters(args)
        for feat, required in step_conv_solv.undefined():
            step_conv_solv.use(self._get(feat | SOP.ForIncr, required))
        self.use(step_conv_solv)
        return step_conv_solv

    def _get_step_solver(self):
        logger.debug("+++ get StepSolver")
        step_solver = self.get_feature(SOP.StepSolver, optional=True)
        if not step_solver:
            step_solver = StepSolver()
        # FIXME replace by 'use(*, Keywords)'
        if not step_solver.param:
            step_solver.setParameters(self.get_feature(SOP.Keywords))
        for feat, required in step_solver.undefined():
            step_solver.use(self._get(feat | SOP.ForStep, required))
        self.use(step_solver)
        return step_solver

    def _get(self, option, required):
        logger.debug(f"--- requesting for {SOP.name(option)}")
        if option & SOP.PhysicalProblem:
            return self.phys_pb
        if option & SOP.PhysicalState:
            return self.phys_state
        if option & SOP.Storage:
            return self._get_storage()
        if option & SOP.LinearSolver:
            return self._get_linear_solver()
        if option & SOP.Contact:
            return self._get_contact_manager()
        if option & (SOP.ConvergenceManager | SOP.ForIncr) == option:
            return self._get_incr_conv()
        if option & SOP.IncrementalSolver:
            return self._get_incremental_solver()
        if option & (SOP.ConvergenceManager | SOP.ForStep) == option:
            return self._get_step_conv()
        if option & SOP.ConvergenceCriteria:
            return self._get_step_conv_solver()
        if option & SOP.StepSolver:
            return self._get_step_solver()
        if required:
            raise NotImplementedError(f"unsupported feature id: {option}")
