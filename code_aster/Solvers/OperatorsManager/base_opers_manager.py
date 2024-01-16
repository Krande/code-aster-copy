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

from ..Basics import SolverFeature
from ..Basics import SolverOptions as SOP
from ...Objects import DiscreteComputation
from ...Utilities import profile


class BaseOperatorsManager(SolverFeature):
    """Solve an iteration."""

    provide = SOP.OperatorsManager

    required_features = [SOP.PhysicalProblem, SOP.PhysicalState, SOP.LinearSolver]

    optional_features = [SOP.Contact]

    def __init__(self):
        super().__init__()

    def initialize(self):
        """Initializes the operator manager."""
        raise NotImplementedError

    def finalize(self):
        """Finalizes the operator manager."""
        raise NotImplementedError

    def executeIteration(self, iter_idx):
        """Should Newton iteration iter_idx be performed

        Arguments:
            iter_idx (int): Newton iteration number.

        Returns:
            bool: whether Newton's iteration should be excuted or
            not, even if the solver has converged
        """
        return False

    @profile
    @SolverFeature.check_once
    def getResidual(self, scaling=1.0):
        """Compute R(u, Lagr) = - (Rint(u, Lagr) + Rcont(u, Lagr) - Rext(u, Lagr)).

        This is not the true residual but the opposite.

        Arguments:
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0).

        Returns:
            tuple (Residuals(), FieldOnCellsReal, FieldOnCellsReal):
            Tuple with residuals, internal state variables (VARI_ELGA),
            Cauchy stress tensor (SIEF_ELGA).
        """

        disc_comp = DiscreteComputation(self.phys_pb)
        contact_manager = self.get_feature(SOP.Contact, optional=True)

        return disc_comp.getResidual(
            self.phys_state, contact_manager=contact_manager, scaling=scaling
        )

    @profile
    @SolverFeature.check_once
    def getStiffnessJacobian(self, matrix_type):
        """Compute K(u) = d(Rint(u) - Rext(u)) / du

        Arguments:
            matrix_type (str): type of matrix used.

        Returns:
            AssemblyMatrixDisplacementReal: Jacobian matrix.
        """
        disc_comp = DiscreteComputation(self.phys_pb)
        contact_manager = self.get_feature(SOP.Contact, optional=True)

        return disc_comp.getTangentMatrix(
            self.phys_state, matrix_type=matrix_type, contact_manager=contact_manager, assemble=True
        )
