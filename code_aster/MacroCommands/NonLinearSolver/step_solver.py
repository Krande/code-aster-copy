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

from libaster import deleteTemporaryObjects

from ...Supervis import ConvergenceError
from ...Utilities import no_new_attributes, profile
from .convergence_manager import ConvergenceManager
from .incremental_solver import IncrementalSolver
from .logging_manager import LoggingManager


class StepSolver:
    """Solves a step, loops on iterations.

    Arguments:
        step_rank (int): Current step to solve.
    """

    current_iter = phys_state = prediction = None
    phys_pb = linear_solver = None
    epsilon = max_iter = epsilon_type = None
    matr_update_iter = matr_update_step = current_matrix = step_rank = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, step_rank):
        self.current_iter = 0
        self.step_rank = step_rank

    def setParameters(self, epsilon_maxi, epsilon_rela, max_iter):
        """Define convergency parameters.

        `epsilon_maxi` is only used if `epsilon_real` is *None*.

        Arguments:
            epsilon_maxi (float): Absolute criteria.
            epsilon_rela (float): Relative criteria.
            max_iter (int): Maximum number of iterations.
        """
        if epsilon_rela is None:
            self.epsilon = epsilon_maxi
            self.epsilon_type = "RESI_GLOB_MAXI"
        else:
            self.epsilon = epsilon_rela
            self.epsilon_type = "RESI_GLOB_RELA"
        self.max_iter = max_iter

    def setPhysicalProblem(self, phys_pb):
        """Assign the physical problem.

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
        """
        self.phys_pb = phys_pb

    def setPhysicalState(self, phys_state):
        """Assign the physical state.

        Arguments:
            phys_state (PhysicalState): Physical state
        """
        self.phys_state = phys_state
        self.phys_state.displ_incr = self.phys_state.createDisplacement(
            self.phys_pb, 0.0)

    def getPhysicalState(self):
        """Get the physical state.

        Returns:
            PhysicalState: Currently used Physical state
        """
        return self.phys_state

    def setLinearSolver(self, linear_solver):
        """Set the linear solver to be used.

        Arguments:
            linear_solver (LinearSolver): Linear solver object.
        """
        self.linear_solver = linear_solver

    def updatePhysicalState(self, displ_incr, variP, sigma, convManager):
        """Update the physical state.

        Arguments:
            displ_incr (FieldOnNodes): Displacement increment.
            variP (FieldOnCells): Internal state variables.
            sigma (FieldOnCells): Stress field.
            convManager (ConvergenceManager): Object that manages the
                convergency criteria.
        """
        self.phys_state.displ_incr += displ_incr

        if convManager.hasConverged():
            self.phys_state.variP = variP
            self.phys_state.stress = sigma

    def setPrediction(self, prediction):
        """Select type of prediction.

        Arguments
            prediction (str): predicition used in "ELASTIQUE" or "TANGENTE"
        """
        assert prediction in (
            "ELASTIQUE", "TANGENTE"), f"unsupported value: {prediction}"
        self.prediction = prediction

    def setUpdateParameters(self, REAC_INCR, REAC_ITER):
        """Set REAC_INCR and REAC_ITER parameters.

        Arguments:
            REAC_INCR (int): Number of steps between updating.
            REAC_ITER (int): Number of iterations between updating.
        """
        self.matr_update_step = REAC_INCR
        self.matr_update_iter = REAC_ITER

    def createConvergenceManager(self):
        """Return an object that holds convergence criteria.

        Returns:
            ConvergenceManager: object that manages the convergency criteria.
        """
        return ConvergenceManager(self.epsilon, self.epsilon_type, self.phys_pb, self.phys_state)

    def createIncrementalSolver(self):
        """Return a solver for the next iteration.

        Returns:
            IncrementalSolver: object to solve the next iteration.
        """
        # current_state == U, Delta_U, P, (Sigm)
        return IncrementalSolver()

    def createLoggingManager(self):
        """Return a logging manager

        Returns:
            LoggingManager: object for logging
        """
        logManager = LoggingManager()
        logManager.addConvTableColumn("NEWTON")
        logManager.addConvTableColumn("RESIDU RELATIF RESI_GLOB_RELA")
        logManager.addConvTableColumn("RESIDU ABSOLU RESI_GLOB_MAXI")
        logManager.addConvTableColumn("OPTION ASSEMBLAGE")

        return logManager

    def hasFinished(self):
        """Tell if there are iterations to be computed.

        Returns:
            bool: *True* if there is no iteration to be computed, *False* otherwise.
        """
        return self.current_iter > self.max_iter

    def _setMatrixType(self):
        """Set matrix type.

        Returns:
            str: Type of matrix to be computed.
        """
        if self.current_iter == 0:
            matrix_type = "PRED_" + self.prediction
            if self.step_rank % self.matr_update_step == 0:
                # make unavailable the current predicted matrix
                self.current_matrix = None
        else:
            matrix_type = "TANGENTE"
            if self.current_iter % self.matr_update_iter == 0:
                # make unavailable the current tangent matrix
                self.current_matrix = None
        return matrix_type

    @profile
    def solve(self):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        convManager = self.createConvergenceManager()
        logManager = self.createLoggingManager()
        logManager.printIntro(self.phys_state.time + self.phys_state.time_step, 1)
        logManager.printConvTableEntries()

        while not self.hasFinished() and not convManager.hasConverged():
            iteration = self.createIncrementalSolver()
            iteration.setConvergenceCriteria(convManager)
            iteration.setPhysicalProblem(self.phys_pb)
            iteration.setPhysicalState(self.phys_state)
            iteration.setLinearSolver(self.linear_solver)

            # select type of matrix
            matrix_type = self._setMatrixType()
            displ_incr, variP, sigma, self.current_matrix = iteration.solve(
                matrix_type, self.current_matrix
            )
            self.updatePhysicalState(displ_incr, variP, sigma, convManager)

            logManager.printConvTableRow([self.current_iter, convManager.residual["RESI_GLOB_RELA"],
                                          convManager.residual["RESI_GLOB_MAXI"], matrix_type])

            self.current_iter += 1

        if not convManager.hasConverged():
            raise ConvergenceError("MECANONLINE9_12")

        # Je pense que l'on garde trop les matrices assemblées et les factorisées.
        # Il faudrait la supprimer si on ne compte pas la garder
        if self.matr_update_step == 1 or self.step_rank % self.matr_update_step != 0:
            self.current_matrix = None

        deleteTemporaryObjects()

        logManager.printConvTableEnd()
