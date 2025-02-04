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

from libaster import deleteTemporaryObjects

from .base_step_solver import BaseStepSolver
from ..Basics import SolverOptions as SOP
from ...Utilities import no_new_attributes, profile
from ..OperatorsManager import MecaStatOperatorsManager
from ..Basics import ProblemType


class MecaStatStepSolver(BaseStepSolver):
    """Solves a step, loops on iterations."""

    problem_type = ProblemType.MecaStat

    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def create(cls, param=None):
        """Setup a solver for the given problem.

        Arguments:
            param (dict) : user keywords.

        Returns:
            *StepSolver*: A relevant *StepSolver* object.
        """
        return cls()

    def __init__(self):
        super(MecaStatStepSolver, self).__init__()

    def initialize(self):
        """Initialization."""
        super().initialize()
        self.phys_state.primal_step = self.phys_state.createPrimal(self.phys_pb, 0.0)

    # @profile
    def solve(self):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        logManager = self.createLoggingManager()
        logManager.printConvTableEntries()

        self.conv_criteria.setLoggingManager(logManager)
        self.conv_criteria.initialize()

        # Solve current iteration
        self.current_matrix = self.conv_criteria.solve(self.current_matrix)

        deleteTemporaryObjects()

        logManager.printConvTableEnd()

    def setup(self):
        """set up the step solver."""
        opers_manager = self.get_feature(SOP.OperatorsManager, optional=True)
        if not opers_manager:
            opers_manager = MecaStatOperatorsManager()
        for feat, required in opers_manager.undefined():
            feat_obj = self.get_feature(feat, optional=(not required))
            opers_manager.use(feat_obj)
        self.conv_criteria.use(opers_manager)
        if hasattr(self.conv_criteria, "local_solver"):
            if self.conv_criteria.local_solver:
                self.conv_criteria.local_solver.use(opers_manager)
        self.use(opers_manager)
