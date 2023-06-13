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
from .geometric_solver import GeometricSolver
from .incremental_solver import IncrementalSolver
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
        SOP.IncrementalSolver,
        SOP.ConvergenceCriteria,
        SOP.ConvergenceManager,
        SOP.ForStep,
        SOP.ForIncr,
    ]
    optional_features = [SOP.Contact, SOP.PostStepHook]

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

    def _get_contact_manager(self):
        logger.debug("+++ get ContactManager")
        return self.get_feature(SOP.Contact, optional=True)

    def _get_incr_conv(self):
        logger.debug("+++ get ConvergenceManager ForIncr")
        converg = self.get_feature(SOP.ConvergenceManager | SOP.ForIncr, optional=True)
        if not converg:
            args = self.get_feature(SOP.Keywords)
            converg = ConvergenceManager()
            for crit in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI", "ITER_GLOB_MAXI"):
                value = args["CONVERGENCE"].get(crit)
                if value is not None:
                    converg.setdefault(crit, value)
            if args.get("CONTACT"):
                if args["CONTACT"].get("ALGO_RESO_GEOM") == "NEWTON":
                    converg.setdefault("RESI_GEOM", args["CONTACT"].get("RESI_GEOM"))
            if not converg.hasResidual():
                converg.setdefault("RESI_GLOB_RELA", 1.0e-6)
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
            nbIterMini = 0
            nbIterMaxi = 1
            if args.get("CONTACT"):
                contact = args["CONTACT"]
                value = 1.0e150
                value = contact.get("RESI_GEOM", value)
                step_crit.setdefault("RESI_GEOM", value)
                if contact.get("REAC_GEOM") == "AUTOMATIQUE":
                    nbIterMini = 0
                    nbIterMaxi = contact["ITER_GEOM_MAXI"]
                elif contact.get("REAC_GEOM") == "CONTROLE":
                    nbIterMaxi = contact["NB_ITER_GEOM"]
                    nbIterMini = nbIterMaxi - 1
            step_crit.setdefault("ITER_GEOM", nbIterMaxi).minValue = nbIterMini
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
