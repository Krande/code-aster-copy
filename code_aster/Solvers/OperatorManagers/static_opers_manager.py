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

from .base_opers_manager import BaseOperatorsManager
from ...Utilities import no_new_attributes


class StaticOperatorsManager(BaseOperatorsManager):
    """Solve an iteration."""

    _first_jacobian = _lagr_scaling = None
    _temp_stress = _temp_internVar = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()
        self._first_jacobian = self._lagr_scaling = None
        self._temp_stress = self._temp_internVar = None

    def initialize(self):
        """Initializes the operator manager."""
        self._first_jacobian = self._lagr_scaling = None
        self._temp_stress = self._temp_internVar = None

    def finalize(self):
        """Finalizes the operator manager."""
        self.phys_state.stress = self._temp_stress
        self.phys_state.internVar = self._temp_internVar

    @property
    def first_jacobian(self):
        """Returns the first computed Jacobian"""
        assert self._first_jacobian is not None
        return self._first_jacobian

    def getLagrangeScaling(self, matrix_type):
        """Returns Lagrange scaling.

        Arguments:
            matrix_type (str): type of matrix used.

        """

        if self._lagr_scaling is None:
            self._first_jacobian = self.getJacobian(matrix_type)

        return self._lagr_scaling

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
        resi_state, internVar, stress = super().getResidual(scaling=scaling)
        self._temp_stress = stress
        self._temp_internVar = internVar
        return resi_state

    def getJacobian(self, matrix_type):
        """Compute K(u) = d(Rint(u) - Rext(u)) / du

        Arguments:
            matrix_type (str): type of matrix used.

        Returns:
            AssemblyMatrixDisplacementReal: Jacobian matrix.
        """
        if self._first_jacobian is not None:
            jacobian = self._first_jacobian
            self._first_jacobian = None
        else:
            jacobian = super().getStiffnessJacobian(matrix_type)
            self._lagr_scaling = jacobian.getLagrangeScaling()
        return jacobian
