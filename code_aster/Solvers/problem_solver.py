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

from libaster import deleteTemporaryObjects, resetFortranLoggingLevel, setFortranLoggingLevel

from ..Objects import LinearSolver
from ..Utilities import DEBUG, INFO, WARNING, ExecutionParameter, Options, logger, no_new_attributes
from .convergence_manager import ConvergenceManager
from .incremental_solver import IncrementalSolver
from .line_search import LineSearch
from .physical_state import PhysicalState
from .snes_solver import SNESSolver
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP
from .step_solver import StepSolver
from .storage_manager import StorageManager


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
        SOP.SnesSolver,
        SOP.ConvergenceManager,
    ]
    optional_features = [
        SOP.Contact,
        SOP.PostStepHook,
        SOP.LineSearch,
        SOP.IncrementalSolver,
        SOP.ConvergenceCriteria,
    ]

    _main = _result = None
    _phys_state = None
    _verb = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, main, result) -> None:
        super().__init__()
        self._main = main
        self._result = result
        self._phys_state = PhysicalState()
        self._verb = logger.getEffectiveLevel(), ExecutionParameter().option & Options.ShowSyntax

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
        args = _F(args)
        self.use(args, SOP.Keywords)

    def initialize(self):
        """Initialize default objects when required."""
        args = self.get_feature(SOP.Keywords, optional=True)
        self._setLoggingLevel(args.get("INFO", 1))
        # required to build other default objects
        self._main.use(self.phys_pb)
        self._main.use(self._phys_state)
        self._main.use(args, SOP.Keywords)
        self._main.use(self.get_feature(SOP.TimeStepper))
        self._main.use(self._get_storage())
        self._main.use(self._get_step_solver())
        self._main.use(self._get_snes_solver())
        self._main.use(self.get_features(SOP.PostStepHook), provide=SOP.PostStepHook)
        self.check_features()

    def run(self):
        """Solve the problem.

        Returns:
            *misc*: result object.
        """
        try:
            self.initialize()
            self._main.run()
        finally:
            self._resetLoggingLevel()
            deleteTemporaryObjects()
        return self._result

    def _setLoggingLevel(self, level):
        """Set logging level.

        Arguments:
            level (int): verbosity level (meaning INFO keyword).
        """
        info = {0: WARNING, 1: INFO, 2: DEBUG, 3: DEBUG, 4: DEBUG}
        if level is None:
            level = 1
        setFortranLoggingLevel(level)
        logger.setLevel(info[level])
        # Disable printing of python command
        if level < 3:
            ExecutionParameter().disable(Options.ShowSyntax)

    def _resetLoggingLevel(self):
        """Reset logging level."""
        level, show = self._verb
        resetFortranLoggingLevel()
        logger.setLevel(level)
        if show:
            ExecutionParameter().enable(Options.ShowSyntax)

    # methods that build and return objects (allowing dependencies)
    def _get_storage(self):
        logger.debug("+++ get StorageManager")
        store = self.get_feature(SOP.Storage, optional=True)
        if not store:
            args = self.get_feature(SOP.Keywords)
            store = StorageManager(self._result, args.get("ARCHIVAGE"))
        self.use(store)
        return store

    def _get_linear_solver(self):
        logger.debug("+++ get LinearSolver")
        if not self.has_feature(SOP.LinearSolver):
            args = self.get_feature(SOP.Keywords)
            self.use(LinearSolver.factory("STAT_NON_LINE", args["SOLVEUR"]))
        return self.get_feature(SOP.LinearSolver)

    def _get_line_search(self):
        logger.debug("+++ get LineSearch")
        if not self.has_feature(SOP.LineSearch):
            args = self.get_feature(SOP.Keywords)
            line = LineSearch(args.get("RECH_LINEAIRE"))
        for feat, required in line.undefined():
            line.use(self._get(feat, required))
        self.use(line)
        return self.get_feature(SOP.LineSearch)

    def _get_contact_manager(self):
        logger.debug("+++ get ContactManager")
        return self.get_feature(SOP.Contact, optional=True)

    def _get_conv_manager(self):
        logger.debug("+++ get ConvergenceManager")
        converg = self.get_feature(SOP.ConvergenceManager, optional=True)
        if not converg:
            args = self.get_feature(SOP.Keywords)
            converg = ConvergenceManager()
            for crit in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI", "ITER_GLOB_MAXI"):
                value = args["CONVERGENCE"].get(crit)
                if value is not None:
                    converg.setdefault(crit, value)
            if not converg.hasResidual():
                converg.setdefault("RESI_GLOB_RELA", 1.0e-6)
            if args.get("CONTACT"):
                converg.setdefault("RESI_GEOM", args["CONTACT"].get("RESI_GEOM"))
        for feat, required in converg.undefined():
            converg.use(self._get(feat, required))
        self.use(converg)
        return converg

    def _get_incremental_solver(self):
        logger.debug("+++ get IncrementalSolver")
        incr_solver = self.get_feature(SOP.IncrementalSolver, optional=True)
        if not incr_solver:
            incr_solver = IncrementalSolver()
        for feat, required in incr_solver.undefined():
            incr_solver.use(self._get(feat, required))
        self.use(incr_solver)
        return incr_solver

    def _get_step_conv_solver(self):
        logger.debug("+++ get ConvergenceCriteria")
        step_conv_solv = self.get_feature(SOP.ConvergenceCriteria, optional=True)
        args = self.get_feature(SOP.Keywords)
        if not step_conv_solv:
            step_conv_solv = self._get(SOP.IncrementalSolver, True)
        for feat, required in step_conv_solv.undefined():
            step_conv_solv.use(self._get(feat, required))
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
            step_solver.use(self._get(feat, required))
        self.use(step_solver)
        return step_solver

    def _get_snes_solver(self):
        logger.debug("+++ get SnesSolver")
        snes_solver = self.get_feature(SOP.SnesSolver, optional=True)
        if not snes_solver:
            snes_solver = SNESSolver()
        # FIXME replace by 'use(*, Keywords)'
        if not snes_solver.param:
            snes_solver.setParameters(self.get_feature(SOP.Keywords))
        for feat, required in snes_solver.undefined():
            snes_solver.use(self._get(feat, required))
        self.use(snes_solver)
        return snes_solver

    def _get(self, option, required):
        if option & SOP.PhysicalProblem:
            return self.phys_pb
        if option & SOP.PhysicalState:
            return self.phys_state
        if option & SOP.Storage:
            return self._get_storage()
        if option & SOP.LinearSolver:
            return self._get_linear_solver()
        if option & SOP.LineSearch:
            return self._get_line_search()
        if option & SOP.Contact:
            return self._get_contact_manager()
        if option & SOP.ConvergenceManager:
            return self._get_conv_manager()
        if option & SOP.IncrementalSolver:
            return self._get_incremental_solver()
        if option & SOP.ConvergenceCriteria:
            return self._get_step_conv_solver()
        if option & SOP.StepSolver:
            return self._get_step_solver()
        if option & SOP.SnesSolver:
            return self._get_snes_solver()
        if required:
            raise NotImplementedError(f"unsupported feature id: {option}")
