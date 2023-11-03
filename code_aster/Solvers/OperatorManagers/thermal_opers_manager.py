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
from ...Objects import DiscreteComputation
from ..residual import Residuals

class ThermalOperatorsManager(BaseOperatorsManager):
    """Solve an iteration."""

    _first_jacobian = _lagr_scaling = None
    _temp_stress = _temp_internVar = _theta = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, theta=None):
        super().__init__()
        self._first_jacobian = self._lagr_scaling = None
        self._temp_stress = self._temp_internVar = None
        self._theta = theta

    def isStationary(self):
        """Checks whether it is a stationary or transient case"""
        return self._theta is None

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
        """Computes the residual."""
        if self.isStationary():
            result = self._getResidualStat(scaling)
        else:
            result = self._getResidualTrans(scaling)
        return result

    def _getResidualStat(self, scaling=1.0):
        """Computes the residual."""
        resi_state, internVar, stress = super().getResidual(scaling=scaling)
        self._temp_stress = stress
        self._temp_internVar = internVar
        return resi_state

    def _getResidualTrans(self, scaling=1.0):
        """Computes the residual."""
        raise NotImplementedError

    def getJacobian(self, matrix_type):
        """Compute the jacobian matrix"""
        if self._first_jacobian is not None:
            jacobian = self._first_jacobian
            self._first_jacobian = None
        else:
            if self.isStationary():
                jacobian = self._getJacobianStat(matrix_type)
            else:
                jacobian = self._getJacobianTrans(matrix_type)
            self._lagr_scaling = jacobian.getLagrangeScaling()
        return jacobian

    def _getJacobianStat(self, matrix_type):
        """Compute the jacobian matrix"""
        return super().getStiffnessJacobian(matrix_type)

    def _getJacobianTrans(self, matrix_type):
        """Compute the jacobian matrix"""
        raise NotImplementedError
