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
from math import sqrt


class OnSubStepIntegrator(BaseIntegrator):
    """
    Integrator that works on a sub-step.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._intrm_state = None

    def setIntermediateState(self, state):
        """Define the state at the beginning of the iteration.

        Arguments:
            state (PhysicalState): State at the beginning of the sub-step.
        """
        self._intrm_state = state.duplicate()

    @property
    def Um(self):
        """Intermediate primal unknowns."""
        return self._intrm_state.current.U

    @Um.setter
    def Um(self, field):
        """Set intermediate primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self._intrm_state.current.U = field

    @property
    def dUm(self):
        """Intermediate derivative of primal unknowns."""
        return self._intrm_state.current.dU

    @dUm.setter
    def dUm(self, field):
        """Set derivative of intermediate primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self._intrm_state.current.dU = field

    @property
    def d2Um(self):
        """Intermediate second derivative of primal unknowns."""
        return self._intrm_state.current.d2U

    @d2Um.setter
    def d2Um(self, field):
        """Set second derivative of intermediate primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self._intrm_state.current.d2U = field


class BDF2Integrator(OnSubStepIntegrator):
    """
    Implementation of Backward Differentiation Formula scheme of order 2
    """

    integration_type = IntegrationType.Implicit
    integrator_name = IntegratorName.Bdf2

    _gamma = _G0 = _dFc = None

    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def create(cls, schema):
        """Setup a time integrator for the given problem.

        Arguments:
            schema (dict) : *SCHEMA_TEMPS* keyword.

        Returns:
            *BDF2Integrator*: A *BDF2Integrator* object.
        """
        return cls()

    def __init__(self, gamma=2 - sqrt(2)):
        super().__init__()
        self._gamma = gamma
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
            g3 = (1 - self._gamma) / (2 - self._gamma)
            jacobian = self._mass / g3**2 / self.dt**2 + C / g3 / self.dt + K
            self._lagr_scaling = jacobian.getLagrangeScaling()
        return jacobian

    def getResidual(self, scaling=1.0):
        """Compute the residue vector."""
        g1 = 1 / self._gamma / (2 - self._gamma)
        g2 = -((1 - self._gamma) ** 2) / self._gamma / (2 - self._gamma)
        g3 = (1 - self._gamma) / (2 - self._gamma)
        G0 = self.U - g1 * self.Um - g2 * self.U0 - g3 * self.dt * self.dU
        self._G0 = G0
        G1 = self.dU - g1 * self.dUm - g2 * self.dU0 - g3 * self.dt * self.d2U
        residue = (self._mass * G1) / g3 / self.dt + (
            self._mass / g3**2 / self.dt**2 - self._dFc / g3 / self.dt
        ) * G0
        return residue

    def updateVariables(self, deltaU):
        """Update the physical state."""
        g3 = (1 - self._gamma) / (2 - self._gamma)
        # self.U += deltaU
        self.dU += (deltaU + self._G0) / g3 / self.dt
        self.computeAcceleration()
