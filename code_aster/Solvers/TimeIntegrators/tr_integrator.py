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

from .base_integrator import BaseIntegrator, IntegrationType, IntegratorName
from ...Utilities import no_new_attributes


class TRIntegrator(BaseIntegrator):
    """Implementation of Trapezoidal rule scheme"""

    integration_type = IntegrationType.Implicit
    integrator_name = IntegratorName.Tr

    _G0 = _dFc = None

    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def create(cls, schema):
        """Setup a time integrator for the given problem.

        Arguments:
            schema (dict) : *SCHEMA_TEMPS* keyword.

        Returns:
            *TRIntegrator*: A *TRIntegrator* object.
        """
        return cls()

    def __init__(self):
        super().__init__()
        self._G0 = self._dFc = None

    def initializeStep(self):
        """Define the step parameters."""
        # A supprimer car redondante
        self.phys_state = self._init_state.duplicate()
        self.computeAcceleration()

    def getJacobian(self, matrix_type):
        """Compute the jacobian matrix.

        Arguments:
            matrix_type (str): type of matrix used.

        Returns:
            AssemblyMatrixDisplacementReal: Jacobian matrix.
        """
        if self._first_jacobian is not None:
            jacobian = self._first_jacobian
            self._first_jacobian = None
        else:
            K, C = self.getStiffAndDamp(self.t0, self.dt, self.U, self.dU, self.d2U, matrix_type)
            self._dFc = C
            jacobian = self._mass * 4.0 / self.dt**2 + 2.0 * C / self.dt + K
            self._lagr_scaling = jacobian.getLagrangeScaling()
        return jacobian

    def getResidual(self, scaling=1.0):
        """Compute the residue vector."""
        self._G0 = self.U - self.U0 - 0.5 * self.dt * (self.dU + self.dU0)
        G1 = self.dU - self.dU0 - 0.5 * self.dt * (self.d2U + self.d2U0)
        residue = (
            self._mass * G1 * 2 / self.dt
            + (self._mass * 4 / self.dt**2 - 2 / self.dt * self._dFc) * self._G0
        )
        return residue

    def updateVariables(self, deltaU):
        """Update the physical state."""
        # self.U += deltaU
        self.dU += 2 * (deltaU + self._G0) / self.dt
        self.computeAcceleration()
