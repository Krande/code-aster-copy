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

from ..TimeIntegrators import TimeScheme, BaseIntegrator
from .meca_dyna_step_solver import MecaDynaStepSolver
from math import sqrt


class MultiStepSolver(MecaDynaStepSolver):
    """Solves a step, loops on iterations."""

    integration_type = TimeScheme.Multiple

    @classmethod
    def create(cls, step_solvers, param):
        """Setup a solver for the given problem.

        Arguments:
            step_solvers (list) : list of step solvers.
            param (dict) : user keywords.

        Returns:
            *StepSolver*: A relevant *StepSolver* object.
        """
        return cls(step_solvers)

    def __init__(self, step_solvers):
        assert len(step_solvers) == 2, len(step_solvers)
        assert isinstance(step_solvers[0]._integrator, BaseIntegrator)
        # assert isinstance(step_solvers[1]._integrator, OnSubStepIntegrator)
        self._step_solvers = step_solvers
        self._coef = 2.0 - sqrt(2)

    def setInitialState(self, initial_state):
        """Define the initial state of the integrator.

        Arguments:
            state (PhysicalState): State at the beginning of the iteration.
        """
        self._step_solvers[0].setInitialState(initial_state)

    def solve(self, t_init, delta_t):
        # supprimer les arguments t_init et delta_t
        step0, step1 = self._step_solvers
        print("++ Solving stage 1")
        state0 = step0.getInitialState()
        step0.solve(t_init, self._coef * delta_t)
        print("++ Solving stage 2")
        step1.setInitialState(state0)
        # step1.setIntermediateState(self.state)
        step1.solve(t_init, delta_t)

    def setup(self):
        """set up the step solver."""
        for step_solver in self._step_solvers:
            for feat, required in step_solver.undefined():
                feat_obj = self.get_feature(feat, optional=(not required))
                step_solver.use(feat_obj)
            step_solver.setup()
