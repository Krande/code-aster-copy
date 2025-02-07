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
from ...Utilities import no_new_attributes, profile
from ..Basics import ContextMixin, EventSource
from ..Basics import SolverOptions as SOP

USE_SCALING = False  # for testing only
PERTURB_JAC = False  # for checking only (sloooow)


class IncrementalSolver(EventSource, ContextMixin):
    """Solve an iteration."""

    _data = None
    _converg = _line_search = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()
        self._data = {}

    # FIXME: to be removed after merging with NewtonSolver
    def share(self, converg, line_search):
        self._converg = converg
        self._line_search = line_search

    def notifyObservers(self, convManager, matrix_type):
        """Notify all observers about the convergence.

        Arguments:
            convManager (ConvergenceManager): Object that holds the criteria values.
            matrix_type (str): Type of matrix used.
        """
        self._data = convManager.getParameters()
        self._data["matrix"] = matrix_type
        self._data["isConverged"] = convManager.isConverged()
        super().notifyObservers()

    def get_state(self):
        """Returns the current residuals to be shared with observers."""
        return SOP.IncrementalSolver, self._data

    # @profile
    def solve(self, matrix_type, matrix=None, force=False):
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
            # import pdb

            # pdb.set_trace()

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
                if self._line_search.activated() and not force:
                    primal_incr = self._line_search.solve(primal_incr, scaling)
        else:
            primal_incr = self.state.createPrimal(self.problem, 0.0)

        # evaluate geometric - convergence
        self._converg.evalGeometricResidual(primal_incr)
        self.notifyObservers(self._converg, matrix_type)

        return primal_incr, jacobian, resi_fields
