# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

from libaster import deleteTemporaryObjects, setFortranLoggingLevel, resetFortranLoggingLevel

from ...Objects import NonLinearResult, PhysicalProblem, LinearSolver
from ...Supervis import ConvergenceError, ExecuteCommand, IntegrationError
from ...Utilities import logger, no_new_attributes, profile
from .physical_state import PhysicalState
from .step_solver import StepSolver
from .storage_manager import StorageManager

"""
Notations to use for mechanics.

We assume that we compute the solution at N time steps (t^1, ..., t^N).
Moreover, we assume that the physical state is given at t^0 (and is entierely known).

The solution is known at (t^0, ..., t^{n-1}) and we search the solution at t^{n}.
At t^{n-1}, this is the previous solution (the variables end with '_prev').
At t^{n}, this is the current solution (the variables end with '_curr').

The [0, I] incremental solutions beetwen t^{n} and t^{n} are the incremental solution
(the variables end with '_incr'). It is typically the solution of Newton or Picard
solver at a given iteration.
The total different between t^{n-1} and t^{n} is the step solution (the variables end
with '_step'). This the sum of incremental solution.

For example for the displacement field and a Newton solver.

- displ^{n-1}: previous displacement (displ_prev)
- displ^{n}: current displacement (displ_curr)
- delta_u = -K^{-1}*R(u): incremental displacement (displ_incr)
- Delta_u = sum_{i=0}^I delta_u^i ( = u^{n} - u^{n-1} at convergence): step
  displacement (displ_incr)
"""


class NonLinearSolver:
    """Main object to solve a non linear problem.
    """
    step_rank = stepper = param = phys_state = None
    phys_pb = None
    linear_solver = None
    storage_manager = None
    current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self.phys_state = PhysicalState()
        self.storage_manager = StorageManager(NonLinearResult())

    def setPhysicalProblem(self, model, material, carael=None):
        """Define the physical problem properties.

        Arguments:
            model (Model): Model object.
            material (MaterialField): Material field.
            carael (ElementaryCharacteristics, optional): Elementary characteristics.
        """
        self.phys_pb = PhysicalProblem(model, material, carael)

    def setKeywords(self, **args):
        """Set parameters from user keywords.

        Arguments:
            args (dict) : user keywords.
        """
        self.param = args

    def setStepper(self, stepper):
        """Set time stepper object.

        Arguments:
            stepper (.stepper.TimeStepper): object to be used.
        """
        self.stepper = stepper

    def setBehaviourProperty(self, keywords):
        """Set keywords for behaviour

        Arguments:
            keywords (dict) : keywords as dict
        """
        self.phys_pb.computeBehaviourProperty(keywords, "NON", "NON", 2)

    def setLinearSolver(self, solver=None, keywords=None):
        """Set linear solver from keywords or directly.

        Arguments:
            solver (LinearSolver): a linear solver object
            keywords (dict) : list of keywords to create linear solver
        """
        if solver:
            self.linear_solver = solver
        elif keywords is not None:
            self.linear_solver = LinearSolver.factory(keywords)
        else:
            raise KeyError(
                "At least one argument is expected: solver or keywords")

    def getStepSolver(self, step_rank):
        """Return a solver for the next step.

        Arguments:
            step_rank (int): index of the step.

        Returns:
            StepSolver: object used to solve the step.
        """
        return StepSolver(step_rank)

    def getResult(self):
        """Get the Result object.

        Returns:
            NonLinearResult: the Result object.
        """
        return self.storage_manager.getResult()

    def hasFinished(self):
        """Tell if there are steps to be computed.

        Returns:
            bool: *True* if there is no step to be computed, *False* otherwise.
        """
        return self.stepper.hasFinished()

    def setLoggingLevel(self, level):
        """Set logging level.

        Arguments:
            level (int): verbosity level.
        """
        if level is not None:
            setFortranLoggingLevel(level)
            # Disable printing of python command
            if level < 3:
                ExecuteCommand.level += 1

    def _resetLoggingLevel(self, level):
        """Reset logging level

        Arguments:
            level (int): verbosity level.
        """
        if level is not None:
            if level < 3:
                ExecuteCommand.level -= 1
        resetFortranLoggingLevel()

    def _storeRank(self, rank, time):
        """Store the current physical state.

        Arguments:
            rank (int): rank index for storing
            time (float): current (pseudo)-time.
        """
        self.storage_manager.storeState(
            rank, time, self.phys_pb, self.phys_state)

    @profile
    def _initializeRun(self):
        """Initialize run"""
        self.phys_state.readInitialState(self.phys_pb, self.param)
        self.step_rank = 0
        self._storeRank(self.step_rank, self.phys_state.time)

    @profile
    def run(self):
        """Solve the problem."""
        self.setLoggingLevel(self._get("INFO"))
        self.phys_pb.computeListOfLoads()
        self.phys_pb.computeDOFNumbering()
        self._initializeRun()

        # Solve nonlinear problem
        while not self.hasFinished():
            time_curr = self.stepper.getNext()
            self.step_rank += 1
            self.phys_state.time_step = time_curr - self.phys_state.time

            solv = self.getStepSolver(self.step_rank)
            solv.setPhysicalProblem(self.phys_pb)
            solv.setPhysicalState(self.phys_state)
            solv.setLinearSolver(self.linear_solver)
            solv.setParameters(epsilon_maxi=self._get("CONVERGENCE", "RESI_GLOB_MAXI"),
                               epsilon_rela=self._get(
                                   "CONVERGENCE", "RESI_GLOB_RELA"),
                               max_iter=self._get("CONVERGENCE", "ITER_GLOB_MAXI"))
            # get Newton parameters
            solv.setPrediction(self._get("NEWTON", "PREDICTION"))
            reac_params = {k: self._get("NEWTON", k, 1)
                           for k in ("REAC_ITER", "REAC_INCR")}
            solv.setUpdateParameters(**reac_params)
            solv.current_matrix = self.current_matrix

            try:
                solv.solve()
                self.phys_state.update(solv.getPhysicalState())
                self._storeRank(self.step_rank, time_curr)
            except (ConvergenceError, IntegrationError) as exc:
                logger.error(exc.message)
                self.stepper.raiseError(exc)
            else:
                self.stepper.completed()
            self.current_matrix = solv.current_matrix

        self._resetLoggingLevel(self._get("INFO"))

        deleteTemporaryObjects()

    def _get(self, keyword, parameter=None, default=None):
        """"Return a keyword value"""
        if parameter is not None:
            assert keyword in self.param
            return self.param.get(keyword).get(parameter, default)

        return self.param.get(keyword, default)
