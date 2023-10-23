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


class RK4Integrator(BaseIntegrator):
    """Implementation of Runge-Kutta order 4 scheme."""

    integration_type = IntegrationType.Explicit
    integrator_name = IntegratorName.Rk4

    _k1 = _k2 = _k3 = _k4 = _first = True

    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def create(cls, schema):
        """Setup a time integrator for the given problem.

        Arguments:
            schema (dict) : *SCHEMA_TEMPS* keyword.

        Returns:
            *RK4Integrator*: A *RK4Integrator* object.
        """
        return cls()

    def __init__(self):
        super().__init__()
        self._k1 = self._k2 = None
        self._k3 = self._k4 = None
        self._first = True

    def integrate(self):
        """Explicit integration."""
        if self._first:
            self._first = False
            self.dU0 *= 0.0
            self.d2U0 *= 0.0

        # A supprimer car redondante
        self.phys_state = self._init_state.duplicate()

        self._k1 = self.d2U0

        self._k2 = self.getFunctional(
            self.t0,
            0.5 * self.dt,
            self.U0 + 0.5 * self.dt * self.dU0,
            self.dU0 + 0.5 * self.dt * self._k1,
            self.d2U0,
        ).resi

        self._k3 = self.getFunctional(
            self.t0,
            0.5 * self.dt,
            self.U0 + 0.5 * self.dt * self.dU0 + 0.25 * self.dt**2 * self._k1,
            self.dU0 + 0.5 * self.dt * self._k2,
            self.d2U0,
        ).resi

        self._k4 = self.getFunctional(
            self.t0,
            self.dt,
            self.U0 + self.dt * self.dU0 + 0.5 * self.dt**2 * self._k2,
            self.dU0 + self.dt * self._k3,
            self.d2U0,
        ).resi

        self.phys_state.U = (
            self.U0 + self.dt * self.dU0 + self.dt**2 * (self._k1 + self._k2 + self._k3) / 6.0
        )

        self.phys_state.dU = (
            self.dU0 + self.dt * (self._k1 + 2.0 * self._k2 + 2.0 * self._k3 + self._k4) / 6.0
        )

        self.phys_state.d2U = self.getFunctional(self.t0, self.dt, self.U, self.dU, self.d2U0).resi

        self._init_state = self.phys_state.duplicate()

        return True
