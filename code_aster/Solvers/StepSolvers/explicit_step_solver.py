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

from ..Basics import SolverOptions as SOP
from .meca_dyna_step_solver import MecaDynaStepSolver
from ..TimeIntegrators import TimeScheme


class ExplicitStepSolver(MecaDynaStepSolver):
    """Solves a step, loops on iterations."""

    integration_type = TimeScheme.Explicit

    @classmethod
    def create(cls, integrator, param):
        """Setup a solver for the given problem.

        Arguments:
            integrator : time integrator
            param (dict) : user keywords.

        Returns:
            *StepSolver*: A relevant *StepSolver* object.
        """
        return cls(integrator)

    def solve(self):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        time, time_step = self.state.time_prev, self.state.time_step

        # mandatory
        criteria = self.get_feature(SOP.ConvergenceCriteria)
        criteria.use(self.get_feature(SOP.OperatorsManager))

        self._integrator.initializeStep(time, time_step)
        self._integrator.integrate()
