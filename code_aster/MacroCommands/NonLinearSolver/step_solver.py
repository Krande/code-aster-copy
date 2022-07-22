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
from .contact_solver import ContactSolver
from .logging_manager import LoggingManager
from ...Objects import DiscreteComputation


class StepSolver:
    """Solves a step, loops on iterations.

    Arguments:
        step_rank (int): Current step to solve.
    """

    current_incr = phys_state = prediction = None
    phys_pb = linear_solver = None
    param = None
    matr_update_incr = matr_update_step = current_matrix = step_rank = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, step_rank):
        self.current_incr = 0
        self.step_rank = step_rank

    def setParameters(self, param):
        """Set parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

        self.prediction = self._get("NEWTON", "PREDICTION")
        assert self.prediction in ("ELASTIQUE", "TANGENTE"), f"unsupported value: {self.prediction}"

        self.matr_update_step = self._get("NEWTON", "REAC_ITER", 1)
        self.matr_update_incr = self._get("NEWTON", "REAC_INCR", 1)

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
        self.phys_state.primal_step = self.phys_state.createPrimal(self.phys_pb, 0.0)

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

    def updatePhysicalState(self, primal_incr, internVar, sigma, convManager):
        """Update the physical state.

        Arguments:
            primal_incr (FieldOnNodes): Displacement increment.
            internVar (FieldOnCells): Internal state variables.
            sigma (FieldOnCells): Stress field.
            convManager (ConvergenceManager): Object that manages the
                convergency criteria.
        """
        self.phys_state.primal_step += primal_incr

        if convManager.hasConverged():
            self.phys_state.internVar = internVar
            self.phys_state.stress = sigma

    def createConvergenceManager(self):
        """Return an object that holds convergence criteria.

        Returns:
            ConvergenceManager: object that manages the convergency criteria.
        """

        convMana = ConvergenceManager(self.phys_pb, self.phys_state)

        criterion = ["RESI_GLOB_RELA", "RESI_GLOB_MAXI"]

        for crit in criterion:
            epsilon = self._get("CONVERGENCE", crit)

            if epsilon is not None:
                convMana.addCriteria(crit, epsilon)

        return convMana

    def createIncrementalSolver(self):
        """Return a solver for the next iteration.

        Returns:
            IncrementalSolver: object to solve the next iteration.
        """

        is_solver = IncrementalSolver()
        is_solver.setPhysicalProblem(self.phys_pb)
        is_solver.setPhysicalState(self.phys_state)
        is_solver.setLinearSolver(self.linear_solver)

        return is_solver

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
        return self.current_incr > self._get("CONVERGENCE", "ITER_GLOB_MAXI")

    def _setMatrixType(self, contact):
        """Set matrix type.

        Returns:
            str: Type of matrix to be computed.
        """
        if self.current_incr == 0:
            matrix_type = "PRED_" + self.prediction
            if self.step_rank % self.matr_update_step == 0:
                # make unavailable the current predicted matrix
                self.current_matrix = None
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")
            if self.current_incr % self.matr_update_incr == 0 or contact.enable:
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
        disc_comp = DiscreteComputation(self.phys_pb)
        contact = ContactSolver(self._get("CONTACT", "DEFINITION"))

        while not self.hasFinished() and not convManager.hasConverged():
            iteration = self.createIncrementalSolver()
            iteration.setConvergenceCriteria(convManager)
            iteration.setContactSolver(contact)

            # pairing for contact
            contact.pairing(self.phys_pb)

            # Select type of matrix
            matrix_type = self._setMatrixType(contact)

            # Solve current iteration
            primal_incr, internVar, sigma, self.current_matrix = iteration.solve(
                matrix_type, self.current_matrix
            )

            # Update physical state
            self.updatePhysicalState(primal_incr, internVar, sigma, convManager)
            contact.update(self.phys_state)

            logManager.printConvTableRow(
                [
                    self.current_incr,
                    convManager.getCriteria("RESI_GLOB_RELA"),
                    convManager.getCriteria("RESI_GLOB_MAXI"),
                    matrix_type,
                ]
            )

            self.current_incr += 1

        if not convManager.hasConverged():
            raise ConvergenceError("MECANONLINE9_12")

        # Je pense que l'on garde trop les matrices assemblées et les factorisées.
        # Il faudrait la supprimer si on ne compte pas la garder
        if self.matr_update_step == 1 or self.step_rank % self.matr_update_step != 0 \
            or contact.enable:
            self.current_matrix = None

        deleteTemporaryObjects()

        logManager.printConvTableEnd()

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        if parameter is not None:
            if keyword in self.param:
                if self.param.get(keyword) is not None:
                    return self.param.get(keyword).get(parameter, default)
                return None
            else:
                return None

        return self.param.get(keyword, default)
