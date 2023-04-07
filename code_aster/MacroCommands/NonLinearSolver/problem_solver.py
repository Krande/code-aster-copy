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
from ...Utilities import no_new_attributes
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
        # required to build other default objects
        self._main.use(self.phys_pb)
        self._main.use(self._phys_state)
        self._main.use(args, SOP.Keywords)
        self._main.use(self.get_feature(SOP.TimeStepper))
        self._main.use(self._get_storage())
        self._main.use(self._get_step_solver())
        # FIXME 'use()' or a global logger
        self._main.setLoggingLevel(args.get("INFO", 1))

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
        if not self.has_feature(SOP.Storage):
            self.use(StorageManager(self._result))
        return self.get_feature(SOP.Storage)

    def _get_linear_solver(self):
        if not self.has_feature(SOP.LinearSolver):
            args = self.get_feature(SOP.Keywords)
            self.use(LinearSolver.factory("STAT_NON_LINE", args["SOLVEUR"]))
        return self.get_feature(SOP.LinearSolver)

    def _get_contact_manager(self):
        return self.get_feature(SOP.Contact, optional=True)

    def _get_incr_conv(self):
        if not self.has_feature(SOP.ConvergenceManager | SOP.ForIncr):
            args = self.get_feature(SOP.Keywords)
            converg = ConvergenceManager()
            converg.use(self.phys_pb)
            converg.use(self.phys_state)
            for crit in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI"):
                epsilon = args["CONVERGENCE"].get(crit)
                if epsilon is not None:
                    converg.addCriteria(crit, epsilon)
            if args.get("CONTACT"):
                if args["CONTACT"].get("ALGO_RESO_GEOM") == "NEWTON":
                    converg.addCriteria("RESI_GEOM", args["CONTACT"].get("RESI_GEOM"))
            self.use(converg, SOP.ForIncr)
        return self.get_feature(SOP.ConvergenceManager | SOP.ForIncr)

    def _get_incremental_solver(self):
        if not self.has_feature(SOP.IncrementalSolver):
            incr_solver = IncrementalSolver()
            incr_solver.use(self.get_feature(SOP.PhysicalProblem))
            incr_solver.use(self.phys_state)
            incr_solver.use(self._get_linear_solver())
            incr_solver.use(self._get_contact_manager())
            incr_solver.use(self._get_incr_conv())
            self.use(incr_solver)
        return self.get_feature(SOP.IncrementalSolver)

    def _get_step_conv(self):
        if not self.has_feature(SOP.ConvergenceManager | SOP.ForStep):
            args = self.get_feature(SOP.Keywords)
            step_crit = ConvergenceManager()
            step_crit.use(self.phys_pb)
            step_crit.use(self.phys_state)
            epsilon = 1.0e150
            epsilon = (args.get("CONTACT") or {}).get("RESI_GEOM", epsilon)
            step_crit.addCriteria("RESI_GEOM", epsilon)
            self.use(step_crit, SOP.ForStep)
        return self.get_feature(SOP.ConvergenceManager | SOP.ForStep)

    def _get_step_conv_solver(self):
        if not self.has_feature(SOP.ConvergenceCriteria):
            args = self.get_feature(SOP.Keywords)
            if args.get("METHODE", "NEWTON") == "NEWTON":
                step_conv_solv = GeometricSolver()
            else:
                step_conv_solv = SNESSolver()
            step_conv_solv.use(self.phys_pb)
            step_conv_solv.use(self.phys_state)
            step_conv_solv.setParameters(args)
            step_conv_solv.use(self._get_contact_manager())
            step_conv_solv.use(self._get_incr_conv())
            step_conv_solv.use(self._get_incremental_solver())
            self.use(step_conv_solv)
        return self.get_feature(SOP.ConvergenceCriteria)

    def _get_step_solver(self):
        if not self.has_feature(SOP.StepSolver):
            args = self.get_feature(SOP.Keywords)
            step_solver = StepSolver()
            step_solver.use(self.phys_pb)
            step_solver.use(self.phys_state)
            step_solver.use(self._get_step_conv_solver())
            step_solver.setParameters(args)
            step_solver.use(self._get_step_conv())
            self.use(step_solver)
        return self.get_feature(SOP.StepSolver)
