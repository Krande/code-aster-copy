# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

from ...Objects import AssemblyMatrixDisplacementReal, DiscreteComputation
from ...Supervis import ConvergenceError
from ...Utilities import PETSc, no_new_attributes, profile
from .contact_manager import ContactManager
from .convergence_manager import ConvergenceManager
from .incremental_solver import IncrementalSolver
from .logging_manager import LoggingManager


class SNESSolver:
    """Solves a step, loops on iterations."""

    current_incr = phys_state = prediction = None
    phys_pb = linear_solver = coordinates = None
    param = logManager = contact_manager = None
    matr_update_incr = current_matrix = step_rank = None
    _snes_solver = _incr_solver = _primal_plus = _primal_incr = None
    _scaling = _options = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self.current_incr = 0

    def setParameters(self, param):
        """Assign parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

        self.prediction = self._get("NEWTON", "PREDICTION")
        assert self.prediction in ("ELASTIQUE", "TANGENTE"), f"unsupported value: {self.prediction}"

        self.matr_update_incr = self._get("NEWTON", "REAC_ITER", 1)

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

    def setContactManager(self, contact_manager):
        """Set contact manager.

        Arguments:
            contact_manager (ContactManager): contact manager
        """
        self.contact_manager = contact_manager

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
        is_solver.setContactManager(self.contact_manager)

        return is_solver

    def _setMatrixType(self):
        """Set matrix type.

        Returns:
            str: Type of matrix to be computed.
        """
        if self.current_incr == 0:
            matrix_type = "PRED_" + self.prediction
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")
        return matrix_type

    def _evalFunction(self, snes, X, F):
        # Get the solution increment from PETSc
        self._primal_incr.fromPetsc(snes.getSolutionUpdate())
        # Increment the solution
        self._primal_incr.applyLagrangeScaling(1 / self._scaling)
        self.phys_state.primal_step += self._primal_incr
        self._primal_plus = self.phys_state.primal + self.phys_state.primal_step
        # Build initial residual
        residual, _, _ = self._incr_solver.computeResidual(self._scaling)
        # Apply Lagrange scaling
        residual.resi.applyLagrangeScaling(1 / self._scaling)
        # Apply DirichletBC into the residual
        disc_comp = DiscreteComputation(self.phys_pb)
        time_curr = self.phys_state.time + self.phys_state.time_step
        diriBCs = disc_comp.getIncrementalDirichletBC(time_curr, self._primal_plus)
        self.current_matrix.applyDirichletBC(diriBCs, residual.resi)
        # Copy to PETSc
        residual.resi.toPetsc().copy(F)

    def _evalJacobian(self, snes, X, J, P):
        if self.current_incr % self.matr_update_incr == 0:
            _matrix = self._incr_solver.computeJacobian(self._setMatrixType())
            self.current_matrix = _matrix
            self._scaling = _matrix.getLagrangeScaling()
            _matrix.toPetsc().copy(P)
        self.current_incr += 1
        if J != P:
            J.assemble()  # matrix-free operator
        return PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    @profile
    def solve(self, current_matrix):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        self.current_matrix = current_matrix
        self._incr_solver = self.createIncrementalSolver()
        self.contact_manager.pairing(self.phys_pb)

        _jac = AssemblyMatrixDisplacementReal(self.phys_pb)
        _disc_comp = DiscreteComputation(self.phys_pb)
        _matr_elem_rigi = _disc_comp.getLinearStiffnessMatrix(time=0.0, with_dual=True)
        _jac.addElementaryMatrix(_matr_elem_rigi)
        _jac.assemble()
        self._scaling = _jac.getLagrangeScaling()
        self.current_matrix = _jac
        p_jac = _jac.toPetsc()

        snes = PETSc.SNES().create()

        # register the function in charge of
        # computing the nonlinear residual
        p_resi = p_jac.getVecRight()
        p_resi.set(0)

        self._primal_plus = self.phys_state.primal_step.duplicate()
        self._primal_incr = self.phys_state.primal_step.duplicate()

        snes.setFunction(self._evalFunction, p_resi)

        snes.setJacobian(self._evalJacobian, p_jac)

        def _monitor(snes, its, fgnorm):
            self.logManager.printConvTableRow([its, " - ", fgnorm, " - ", self._setMatrixType()])

        snes.setMonitor(_monitor)

        OptDB = PETSc.Options()
        if not self._options:
            self.linear_solver.build()
            self._options = self.linear_solver.getPetscOptions()
        OptDB.insertString(self._options)
        snes.setFromOptions()

        # solve the nonlinear problem
        b, x = None, p_resi.duplicate()
        x.set(0)  # zero initial guess
        snes.solve(b, x)

        _, internVar, sigma = self._incr_solver.computeResidual()

        self.phys_state.stress = sigma
        self.phys_state.internVar = internVar

        if snes.getConvergedReason() < 0:
            raise ConvergenceError("MECANONLINE9_7")

        return self.current_matrix

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        if parameter is not None:
            if keyword in self.param and self.param.get(keyword) is not None:
                return self.param.get(keyword).get(parameter, default)
            return default

        return self.param.get(keyword, default)
