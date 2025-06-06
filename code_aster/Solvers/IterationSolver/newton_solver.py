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

import numpy as np

from ...LinearAlgebra import MatrixScaler
from ...Objects import DiscreteComputation
from ...Supervis import ConvergenceError
from ...Utilities import no_new_attributes, profile
from ..Basics import EventId, EventSource
from .convergence_manager import ConvergenceManager
from .iteration_solver import BaseIterationSolver
from .line_search import BaseLineSearch

USE_SCALING = False  # for testing only
PERTURB_JAC = False  # for checking only (sloooow)


class NewtonSolver(BaseIterationSolver, EventSource):
    """Solves a step, loops on iterations."""

    __needs__ = ("problem", "state", "keywords", "oper", "linear_solver", "contact")
    solver_type = BaseIterationSolver.SubType.Newton
    _data = _converg = _line_search = None
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
        return instance

    def __init__(self):
        super().__init__()
        self._data = {}

    def initialize(self):
        """Initialize the object for the next step."""
        super().initialize()
        self._converg.initialize()
        iter_glob = self._converg.setdefault("ITER_GLOB_MAXI")
        iter_glob.minValue = 1

    def update(self, primal_incr, resi_fields=None, callback=None):
        """Update the physical state.

        Arguments:
            primal_incr (FieldOnNodes): Displacement increment.
            resi_fields (dict of FieldOnNodes): Fields of residual values
        """

        self.state.primal_step += primal_incr

        for key, field in resi_fields.items():
            self.state.set(key, field)

        if callback:
            callback(primal_incr)

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
            primal_incr, self.current_matrix, resi_fields = self.solve_iteration(
                matrix_type, self.current_matrix, force
            )

            # Update
            self.update(primal_incr, resi_fields, callback)

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
        # we need the matrix to have scaling factor for Lagrange

        if self.contact:
            self.contact.update(self.state)
            self.contact.pairing()

        if not matrix:
            jacobian = self.oper.getJacobian(matrix_type)
        else:
            jacobian = matrix

        # Get scaling for Lagrange (j'aimerais bien m'en débarasser de là)
        scaling = self.oper.getLagrangeScaling(matrix_type)

        # compute residual
        residuals = self.oper.getResidual(scaling)

        # ---------------------------------------------------------------------
        if PERTURB_JAC:
            neq = self.state.primal_step.size()
            jac = np.zeros((neq, neq))
            primal_save = self.state.primal_step.copy()
            res_0 = np.array(residuals.resi.getValues())
            eps = 1.0e-6
            for i_eq in range(neq):
                if abs(self.state.primal_step[i_eq]) < eps:
                    self.state.primal_step[i_eq] = -eps
                    dd = eps
                else:
                    self.state.primal_step[i_eq] *= 1 - eps
                    dd = -self.state.primal_step[i_eq] + primal_save[i_eq]
                res_p = np.array(self.oper.getResidual(scaling).resi.getValues())
                jac[:, i_eq] = (res_p - res_0)[:] / dd
                # return to initial value
                self.state.primal_step = primal_save.copy()
            residuals = self.oper.getResidual(scaling)
            with np.printoptions(precision=3, suppress=False, linewidth=2000):
                print("perturb_jac=", flush=True)
                for row in range(0, neq):
                    print(jac[row, :], flush=True)
        # ---------------------------------------------------------------------

        # evaluate convergence
        resi_fields = self._converg.evalNormResidual(residuals)

        if not self._converg.isConverged() or force:
            disc_comp = DiscreteComputation(self.problem)

            # Compute Dirichlet BC:=
            diriBCs = disc_comp.getIncrementalDirichletBC(
                self.state.time_curr, self.state.primal_curr
            )

            # Solve linear system
            if USE_SCALING:
                S = MatrixScaler.MatrixScaler()
                S.computeScaling(jacobian, merge_dof=[["DX", "DY"], ["LAGS_C", "LAGS_F1"]])
                S.scaleMatrix(jacobian)
                S.scaleRHS(residuals.resi)
            # ------------------------------------------------------------------------
            if PERTURB_JAC:
                A = jacobian.toNumpy()
                neq = A.shape[0]
                neql = range(0, neq // 2)
                with np.printoptions(precision=3, suppress=False, linewidth=2000):
                    print("residual=", np.asarray(residuals.resi.getValues())[neql], flush=True)
                    print("pert_jac vs jac=", flush=True)
                    for row in neql:
                        print(jac[row, neql], flush=True)
                        print(A[row, neql], flush=True)
                        print("-" * 20)
            # ------------------------------------------------------------------------
            if not jacobian.isFactorized():
                self.linear_solver.factorize(jacobian, raiseException=True)
            primal_incr = self.linear_solver.solve(residuals.resi, diriBCs)
            if USE_SCALING:
                S.unscaleSolution(primal_incr)
            # Use line search
            if not self._converg.isPrediction():
                if self._line_search.isEnabled() and not force:
                    primal_incr = self._line_search.solve(primal_incr, scaling)
        else:
            primal_incr = self.state.createPrimal(self.problem, 0.0)

        # evaluate geometric - convergence
        self._converg.evalGeometricResidual(primal_incr)
        self.notifyObservers(matrix_type)

        return primal_incr, jacobian, resi_fields

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
        return EventId.IterationSolver, self._data
