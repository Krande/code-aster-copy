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
from .contact_manager import ContactManager
from .logging_manager import LoggingManager
from ...Objects import DiscreteComputation


class GeometricSolver:
    """Solves a step, loops on iterations.

    """

    current_incr = phys_state = prediction = None
    phys_pb = linear_solver = coordinates =None
    param = geom = logManager = None
    matr_update_incr = current_matrix = step_rank = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self.current_incr = 0

    def setParameters(self, param):
        """Set parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

        self.prediction = self._get("NEWTON", "PREDICTION")
        assert self.prediction in ("ELASTIQUE", "TANGENTE"), f"unsupported value: {self.prediction}"

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

    def setLoggingManager(self, logManager):
        """Assign the logging manager.

        Arguments:
            logManager (LoggingManager): Logging manager.
        """
        self.logManager = logManager

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

    def setCoordinates(self, coordinates):
        """Set the current coordinates to be used.

        Arguments:
            coordinates (MeshCoordinatesField): current coordinates.
        """
        self.geom = coordinates

    def update(self, primal_incr, internVar, sigma, convManager, contManager):
        """Update the physical state.

        Arguments:
            primal_incr (FieldOnNodes): Displacement increment.
            internVar (FieldOnCells): Internal state variables.
            sigma (FieldOnCells): Stress field.
            convManager (ConvergenceManager): Object that manages the
                convergency criteria.
        """

        convManager.evalGeometricResidual(primal_incr)

        self.phys_state.primal_step += primal_incr

        if convManager.hasConverged():
            self.phys_state.internVar = internVar
            self.phys_state.stress = sigma
        elif self._get("CONTACT", "ALGO_RESO_GEOM") == "NEWTON":
            contManager.update(self.phys_state)
            contManager.pairing(self.phys_pb)

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

        if self._get("CONTACT", "ALGO_RESO_GEOM") == "NEWTON":
            epsilon = self._get("CONTACT", "RESI_GEOM")

            convMana.addCriteria("RESI_GEOM", epsilon)

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

    def hasFinished(self):
        """Tell if there are iterations to be computed.

        Returns:
            bool: *True* if there is no iteration to be computed, *False* otherwise.
        """
        return self.current_incr > self._get("CONVERGENCE", "ITER_GLOB_MAXI")

    def _setMatrixType(self, contManager):
        """Set matrix type.

        Returns:
            str: Type of matrix to be computed.
        """
        if self.current_incr == 0:
            matrix_type = "PRED_" + self.prediction
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")
            if self.current_incr % self.matr_update_incr == 0 or contManager.enable:
                # make unavailable the current tangent matrix
                self.current_matrix = None
        return matrix_type

    @profile
    def solve(self, current_matrix):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        self.current_matrix = current_matrix
        convManager = self.createConvergenceManager()
        contManager = ContactManager(self._get("CONTACT", "DEFINITION"), self.geom)
        contManager.pairing(self.phys_pb)

        while not self.hasFinished() and not convManager.hasConverged():
            iteration = self.createIncrementalSolver()
            iteration.setConvergenceCriteria(convManager)
            iteration.setContactManager(contManager)

            # Select type of matrix
            matrix_type = self._setMatrixType(contManager)

            # Solve current iteration
            primal_incr, internVar, sigma, self.current_matrix = iteration.solve(
                matrix_type, self.current_matrix
            )

            # Update
            self.update(primal_incr, internVar, sigma, convManager, contManager)

            self.logManager.printConvTableRow(
                [
                    self.current_incr,
                    convManager.getCriteria("RESI_GLOB_RELA"),
                    convManager.getCriteria("RESI_GLOB_MAXI"),
                    convManager.getCriteria("RESI_GEOM"),
                    matrix_type,
                ]
            )

            self.current_incr += 1

        if not convManager.hasConverged():
            raise ConvergenceError("MECANONLINE9_7")

        if self.current_incr % self.matr_update_incr == 0 or contManager.enable:
            self.current_matrix = None

        deleteTemporaryObjects()

        return self.current_matrix

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        if parameter is not None:
            if keyword in self.param and self.param.get(keyword) is not None:
                return self.param.get(keyword).get(parameter, default)
            else:
                return default

        return self.param.get(keyword, default)
