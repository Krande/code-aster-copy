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

from ..Basics import ContextMixin, LoggingManager, ProblemTypeMixin
from ..ProblemSolver import IterationsSolver


class BaseStepSolver(ABC, ContextMixin, ProblemTypeMixin):
    """Solves a step, loops on iterations."""

    # required_features = [
    #     SOP.PhysicalProblem,          x
    #     SOP.PhysicalState,            x
    #     SOP.ConvergenceCriteria,      x
    #     SOP.LinearSolver,
    # ]
    _iter_solv = None
    current_matrix = None

    def __init__(self):
        super().__init__()
        self.current_matrix = None

    def initialize(self):
        """Initialization."""
        self._iter_solv = IterationsSolver.factory(self.context)

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

    # @abstractmethod
    # def setup(self):
    #     """set up the step solver."""
