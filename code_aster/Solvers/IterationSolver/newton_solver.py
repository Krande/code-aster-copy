# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

import numpy as np

from ...LinearAlgebra import MatrixScaler
from ...Objects import DiscreteComputation
from ...Supervis import ConvergenceError
from ...Utilities import no_new_attributes, profile
from ..Basics import EventId, EventSource
from .convergence_manager import ConvergenceManager
from .iteration_solver import BaseIterationSolver
from .line_search import BaseLineSearch

#  Debug parameters
USE_SCALING = False  # for testing only


#  Newton solver class
class NewtonSolver(BaseIterationSolver, EventSource):
    """Solves a step, loops on iterations."""

    __needs__ = ("problem", "state", "keywords", "oper", "linear_solver", "contact")
    solver_type = BaseIterationSolver.SubType.Newton
    _eventid = EventId.IterationSolver
    _data = _converg = _line_search = None
    _use_scaling = None
    _context = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def builder(cls, context):
        """Builder of NewtonSolver object.

        Args:
            context (Context): Context of the problem.

        Returns:
            instance: New object.
        """
        instance = super().builder(context)
        instance._converg = ConvergenceManager.builder(context)
        instance._line_search = BaseLineSearch.factory(context)
        instance.add_observer(context.stepper)
        instance.add_observer(context.state)
        return instance

    def __init__(self):
        super().__init__()
        self._data = {}
        # - for debug
        self._use_scaling = USE_SCALING
        self._S = None

    def initialize(self):
        """Initialize the object for the next step."""
        super().initialize()
        self._converg.initialize()
        iter_glob = self._converg.setdefault("ITER_GLOB_MAXI")
        iter_glob.minValue = 1

    def update(self, deltaU, resi_fields=None, callback=None):
        """Update the physical state.

        Arguments:
            deltaU (FieldOnNodes): Displacement increment.
            resi_fields (dict of FieldOnNodes): Fields of residual values
        """

        self.state.deltaU += deltaU

        for key, field in resi_fields.items():
            self.state.set(key, field)

        if callback:
            callback(deltaU)

    def _resetMatrix(self):
        """Reset matrix if needed

        Arguments:
            current_incr (int): index of the current increment
        """
        if (
            self.update_matr_incr > 0 and self.current_incr % self.update_matr_incr == 0
        ) or self.contact:
            # make unavailable the current tangent matrix
            self.current_matrix = None

    @profile
    def solve(self, current_matrix, callback=None):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        self.current_matrix = current_matrix

        iter_glob = self._converg.setdefault("ITER_GLOB_MAXI")

        self.oper.initialize()
        while not self._converg.isFinished():
            iter_glob.value = self.current_incr

            # Select type of matrix
            matrix_type = self.matrix_type
            if self.current_incr > 0:
                self._resetMatrix()

            # Should the iteration be executed even if the solver converged ?
            force = self.oper.shouldExecuteIteration(self.current_incr)

            # Solve current iteration
            deltaU, self.current_matrix, resi_fields = self.solve_iteration(
                matrix_type, self.current_matrix, force
            )

            # Update
            self.update(deltaU, resi_fields, callback)
            self.notifyObservers(matrix_type)

            if self.current_incr > 0:
                self.logManager.printConvTableRow(
                    [
                        self.current_incr - 1,
                        self._converg.get("RESI_GLOB_RELA"),
                        self._converg.get("RESI_GLOB_MAXI"),
                        self._converg.get("RESI_REFE_RELA"),
                        self._converg.get("RESI_GEOM"),
                        matrix_type,
                    ]
                )
            self.current_incr += 1

        if not self._converg.isConverged():
            raise ConvergenceError("MECANONLINE9_7")

        self.oper.finalize()
        self._resetMatrix()
        return self.current_matrix

    @profile
    def solve_iteration(self, matrix_type, matrix=None, force=False):
        """Solve the iteration.

        Arguments:
            matrix_type (str): type of matrix used.
            matrix (AssemblyMatrixDisplacementReal, optional): Stiffness matrix
                to be reused.

        Returns:
            tuple (FieldOnNodesReal, AssemblyMatrixDisplacementReal): Tuple
            with incremental primal, Jacobian matrix used (if computed).
        """
        # -- We need the matrix to have scaling factor for Lagrange

        # Specific operations for contact
        self._update_contact_geometry()
        # Preparation of the linear system
        jacobian = self._get_jacobian_iteration(matrix, matrix_type)
        # Get scaling for Lagrange (j'aimerais bien m'en débarasser de là)
        scaling = self.oper.getLagrangeScaling(matrix_type)
        # Compute residual
        residuals = self._compute_residuals(scaling)
        # Evaluate convergence
        resi_fields = self._converg.evalNormResidual(residuals)
        deltaU = self._compute_primal_incr(residuals, jacobian, force, scaling)
        # evaluate geometric - convergence
        self._converg.evalGeometricResidual(deltaU)
        self.notifyObservers(matrix_type)

        return deltaU, jacobian, resi_fields

    def _compute_primal_incr(self, residuals, jacobian, force, scaling):
        """
        Solve the linear system to obtain the current solution increment
        """
        if not self._converg.isConverged() or force:
            disc_comp = DiscreteComputation(self.problem)

            # Compute Dirichlet BC
            diriBCs = disc_comp.getIncrementalDirichletBC(
                self.state.time_curr, self.state.primal_curr
            )

            with MatrixScaler.matrixScaler(
                jacobian,
                residuals.resi,
                merge_dof=[["DX", "DY"], ["LAGS_C", "LAGS_F1"]],
                scaling_type=self._use_scaling,
                verbose=False,
            ) as matScaler:
                # Solve linear system
                if not jacobian.isFactorized():
                    self.linear_solver.factorize(jacobian, raiseException=True)
                primal_incr = self.linear_solver.solve(residuals.resi, diriBCs)
                # Unscale if MatrixScaler is used
                matScaler.unscaleSolution(primal_incr)
                if not self._converg.isPrediction():
                    if self._line_search.isEnabled() and not force:
                        primal_incr = self._line_search.solve(primal_incr, scaling)
        else:
            deltaU = self.state.createPrimal(self.problem, 0.0)

        return deltaU

    def _get_jacobian_iteration(self, matrix, matrix_type):
        """Obtain the Jacobian matric (recompute or use the one provided)

        Arguments:
            matrix_type (str): type of matrix used.
            matrix (AssemblyMatrixDisplacementReal, optional): Stiffness matrix
                to be reused.

        Returns:
            jacobian (AssemblyMatrixDisplacementReal) : Stiffness matrix to use
            in the current iteration
        """
        if not matrix:
            jacobian = self.oper.getJacobian(matrix_type)
        else:
            jacobian = matrix
        return jacobian

    def _compute_residuals(self, scaling):
        """Computation of the residual"""
        residuals = self.oper.getResidual(scaling)
        return residuals

    def _update_contact_geometry(self):
        """Update the pairing for contact"""
        if self.contact:
            self.contact.update(self.state)
            self.contact.pairing()

    def notifyObservers(self, matrix_type):
        """Notify observers about the convergence.

        Arguments:
            matrix_type (str): Type of matrix used.
        """
        self._data = self._converg.getParameters()
        self._data["matrix"] = matrix_type
        self._data["isConverged"] = self._converg.isConverged()
        super().notifyObservers()

    def get_state(self):
        """Returns the current residuals to be shared with observers."""
        return self._data
