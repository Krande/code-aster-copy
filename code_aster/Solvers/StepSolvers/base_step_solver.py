# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

from abc import ABC, abstractmethod

from ..Basics import (
    ContextMixin,
    LoggingManager,
    ProblemDispatcher,
    ProblemTypeMixin,
    SolverFeature,
)
from ..Basics import SolverOptions as SOP


# FIXME: add ABC after removing SolverFeature
class BaseStepSolver(SolverFeature, ProblemDispatcher, ContextMixin, ProblemTypeMixin):
    """Solves a step, loops on iterations."""

    provide = SOP.StepSolver

    required_features = [
        SOP.PhysicalProblem,
        SOP.PhysicalState,
        SOP.ConvergenceCriteria,
        SOP.LinearSolver,
    ]

    optional_features = [SOP.Contact, SOP.OperatorsManager]

    current_matrix = param = None

    def __init__(self):
        super(BaseStepSolver, self).__init__()
        self.current_matrix = self.param = None

    @property
    def conv_criteria(self):
        """ConvergenceCriteria: Convergence criteria object."""
        return self.get_feature(SOP.ConvergenceCriteria)

    def _createPrivate(self, param):
        self.setParameters(param)

    def setParameters(self, param):
        """Set parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

    def initialize(self):
        """Initialization."""
        self.check_features()

    def createLoggingManager(self):
        """Return a logging manager

        Returns:
            LoggingManager: object for logging
        """
        logManager = LoggingManager()
        logManager.addConvTableColumn("NEWTON")
        logManager.addConvTableColumn("RESIDU RELATIF RESI_GLOB_RELA")
        logManager.addConvTableColumn("RESIDU ABSOLU RESI_GLOB_MAXI")
        logManager.addConvTableColumn("RESIDU GEOMETRIQUE RESI_GEOM")
        logManager.addConvTableColumn("OPTION ASSEMBLAGE")

        return logManager

    @abstractmethod
    def solve(self):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """

    def _get(self, keyword, parameter=None, default=None):
        """Return a keyword value"""
        if parameter is not None:
            if keyword in self.param and self.param.get(keyword) is not None:
                return self.param.get(keyword).get(parameter, default)
            else:
                return default

        return self.param.get(keyword, default)

    @abstractmethod
    def setup(self):
        """set up the step solver."""
