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

from .base_step_solver import BaseStepSolver
from ..TimeIntegrators import TimeScheme, IntegratorType, BaseIntegrator
from ..Basics import ProblemType


class MecaDynaStepSolver(BaseStepSolver):
    """Solves a step, loops on iterations."""

    problem_type = ProblemType.MecaDyna
    integration_type = TimeScheme.Unset

    @classmethod
    def create(cls, param=None):
        """Setup a solver for the given problem.

        Arguments:
            param (dict) : user keywords.

        Returns:
            *StepSolver*: A relevant *StepSolver* object.
        """

        def _setup_one(integrator):
            for klass in cls.__subclasses__():
                if klass.support(integrator):
                    return klass.create(integrator, param)

        schema = param.get("SCHEMA_TEMPS")

        integrator_names = cls._getIntegrators(schema)
        solvers_list = []

        for name in integrator_names:
            integrator = BaseIntegrator.create(name, schema)
            solvers_list.append(_setup_one(integrator))

        if len(solvers_list) > 2:
            raise ValueError("only two steps supported yet")

        if len(solvers_list) > 1:
            for klass in cls.__subclasses__():
                if klass.integration_type == TimeScheme.Multiple:
                    solver = klass.create(solvers_list, param)
        else:
            solver = solvers_list[0]

        return solver

    def __init__(self, integrator):
        super(MecaDynaStepSolver, self).__init__()
        self._integrator = integrator

    def initialize(self):
        """Initialization."""
        super().initialize()
        self.state.primal_step = self.state.createPrimal(self.problem, 0.0)
        self.setInitialState(self.state)

    def setInitialState(self, initial_state):
        """Define the initial state of the integrator.

        Arguments:
            initial_state (PhysicalState): State at the beginning of the iteration.
        """
        return self._integrator.setInitialState(initial_state)

    def getInitialState(self):
        """Return the physical state used at the beginning of the iteration.

        Returns:
            PhysicalState: Initial physical state.
        """
        return self._integrator.getInitialState()

    def solve(self):
        """Solve a step."""
        raise NotImplementedError

    @classmethod
    def support(cls, integrator):
        """Tell if the StepSolver supports this kind of integrator.

        Arguments:
            integrator (Integrator): *Integrator* object.

        Returns:
            bool: *True* if the integrator is supported, *False* otherwise.
        """
        return cls.integration_type == integrator.integration_type

    @classmethod
    def _getIntegrators(cls, schema):
        """Returns the list of integrators to use.

        Arguments:
            schema (dict) : *SCHEMA_TEMPS* keyword.

        Returns:
            List[IntegratorType]: list of integrator names.
        """
        assert schema["SCHEMA"] == "NEWMARK", schema["SCHEMA"]
        integrator_name = IntegratorType.Newmark
        return [integrator_name]

    def setup(self):
        """set up the step solver."""
        for feat, required in self._integrator.undefined():
            feat_obj = self.get_feature(feat, optional=(not required))
            self._integrator.use(feat_obj)
        # self._integrator.setup()
        self.use(self._integrator)
