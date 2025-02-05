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

from ...Objects import DiscreteComputation
from ...Utilities import no_new_attributes, profile
from ..Basics import EventSource, SolverFeature
from ..Basics import SolverOptions as SOP


class IncrementalSolver(SolverFeature, EventSource):
    """Solve an iteration."""

    provide = SOP.IncrementalSolver | SOP.EventSource

    required_features = [SOP.PhysicalProblem, SOP.PhysicalState, SOP.LinearSolver, SOP.LineSearch]

    optional_features = [SOP.Contact, SOP.ConvergenceManager, SOP.OperatorsManager]

    _data = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()
        self._data = {}

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

    @property
    def opers_manager(self):
        """OperatorsManager: Operators manager object."""
        return self.get_feature(SOP.OperatorsManager)

    @profile
    @SolverFeature.check_once
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

        contact_manager = self.get_feature(SOP.Contact, optional=True)

        if contact_manager:
            contact_manager.pairing(self.phys_pb)
            contact_manager.update(self.phys_state)

        if not matrix:
            jacobian = self.opers_manager.getJacobian(matrix_type)
        else:
            jacobian = matrix

        # Get scaling for Lagrange (j'aimerais bien m'en débarasser de là)
        scaling = self.opers_manager.getLagrangeScaling(matrix_type)

        # compute residual
        residuals = self.opers_manager.getResidual(scaling)

        # evaluate convergence
        convManager = self.get_feature(SOP.ConvergenceManager)
        resi_fields = convManager.evalNormResidual(residuals)

        if not convManager.isConverged() or force:
            disc_comp = DiscreteComputation(self.phys_pb)

            # Compute Dirichlet BC:=
            diriBCs = disc_comp.getIncrementalDirichletBC(
                self.phys_state.time_curr, self.phys_state.primal_curr
            )
            # import pdb

            # pdb.set_trace()

            # Solve linear system
            linear_solver = self.get_feature(SOP.LinearSolver)
            if not jacobian.isFactorized():
                linear_solver.factorize(jacobian, raiseException=True)
            primal_incr = linear_solver.solve(residuals.resi, diriBCs)

            # Use line search
            lineSearch = self.get_feature(SOP.LineSearch)
            lineSearch.use(self.opers_manager)

            if not convManager.isPrediction():
                if lineSearch.activated() and not force:
                    primal_incr = lineSearch.solve(primal_incr, scaling)
        else:
            primal_incr = self.phys_state.createPrimal(self.phys_pb, 0.0)

        # evaluate geometric - convergence
        convManager.evalGeometricResidual(primal_incr)
        self.notifyObservers(convManager, matrix_type)

        return primal_incr, jacobian, resi_fields
