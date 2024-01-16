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

from libaster import AssemblyMatrixTemperatureReal

from .base_opers_manager import BaseOperatorsManager
from ...Utilities import no_new_attributes
from ...Objects import DiscreteComputation
from ..Basics import Residuals

class ThermalOperatorsManager(BaseOperatorsManager):
    """Solve an iteration."""

    _first_jacobian = _lagr_scaling = None
    _temp_stress = _temp_internVar = _theta = None
    _resi_prev = _resi_temp = _curr_iter = None

    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, theta=None):
        super().__init__()
        self._first_jacobian = self._lagr_scaling = None
        self._temp_stress = self._temp_internVar = None
        self._theta, self._curr_iter = theta, 0
        self._resi_prev, self._resi_temp = None, []

    def isStationary(self):
        """Checks whether it is a stationary or transient case"""
        return self._theta is None

    def initialize(self):
        """Initializes the operator manager."""
        self._first_jacobian = self._lagr_scaling = None
        self._temp_stress = self._temp_internVar = None
        if self._curr_iter == 0:
            disc_comp = DiscreteComputation(self.phys_pb)
            self._resi_prev, _, _ = super().getResidual(scaling=1.0)
            self._resi_prev.resi_mass = disc_comp.getNonLinearCapacityForces(
                self.phys_state.primal_prev,
                self.phys_state.primal_step,
                self.phys_state.time_step,
                self.phys_state.externVar
            )
            self._curr_iter = self._curr_iter + 1

    def finalize(self):
        """Finalizes the operator manager."""
        self.phys_state.stress = self._temp_stress
        self.phys_state.internVar = self._temp_internVar
        self.phys_state.primal_curr = self.phys_state.primal_step
        self._resi_prev = self._resi_temp[-1]

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

    def getJacobian(self, matrix_type):
        """Compute the jacobian matrix."""
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
        """Computes the jacobian matrix for the stationary case."""
        return super().getStiffnessJacobian(matrix_type)

    def _getJacobianTrans(self, matrix_type):
        """Computes the jacobian matrix for the transient case."""
        disc_comp = DiscreteComputation(self.phys_pb)

        mass_ther = disc_comp.getTangentCapacityMatrix(
            self.phys_state.primal_prev,
            self.phys_state.primal_step,
            self.phys_state.externVar,
        )

        rigi_ther = disc_comp.getTangentConductivityMatrix(
            self.phys_state.primal_prev,
            self.phys_state.primal_step,
            self.phys_state.externVar,
            with_dual=False
        )

        rigi_ther_dual = disc_comp.getDualLinearConductivityMatrix()
        dt, theta = self.phys_state.time_step, self._theta


        rigi_ther_ext = disc_comp.getThermalExchangeMatrix(self.phys_state.time_curr)
        rigi_ther_ext *= -1.0

        rigi_ther_ext.addElementaryTerm(
            disc_comp.getThermalTangentNonLinearNeumannMatrix(
                self.phys_state.primal_curr, self.phys_state.time_curr
            ).getElementaryTerms()
        )
        rigi_ther_ext.addElementaryTerm(
            disc_comp.getThermalTangentNonLinearVolumetricMatrix(
                self.phys_state.primal_curr, self.phys_state.time_curr
            ).getElementaryTerms()
        )
        rigi_ther_ext.build()

        result = AssemblyMatrixTemperatureReal(self.phys_pb)
        result.addElementaryMatrix(rigi_ther, theta)
        result.addElementaryMatrix(rigi_ther_ext, -theta)
        result.addElementaryMatrix(rigi_ther_dual)
        result.addElementaryMatrix(mass_ther, 1.0 / dt)
        result.assemble()

        return result

    def getResidual(self, scaling=1.0):
        """Computes the residual."""
        if self.isStationary():
            result = self._getResidualStat(scaling)
        else:
            result = self._getResidualTrans(scaling)
        return result

    def _getResidualStat(self, scaling=1.0):
        """Computes the residual for the stationary case."""
        resi_state, internVar, stress = super().getResidual(scaling=scaling)
        self._temp_stress = stress
        self._temp_internVar = internVar
        return resi_state

    def _getResidualTrans(self, scaling=1.0):
        """Computes the residual for the transient case."""
        dt, theta = self.phys_state.time_step, self._theta
        disc_comp = DiscreteComputation(self.phys_pb)

        temp_prev = self.phys_state.primal_prev
        temp_time = self.phys_state.time_prev

        self.phys_state.primal_prev = self.phys_state.primal_curr
        self.phys_state.time_prev = self.phys_state.time_curr

        resi_curr, internVar, stress = super().getResidual(scaling=scaling)

        self.phys_state.primal_prev = temp_prev
        self.phys_state.time_prev = temp_time

        resi_mass = disc_comp.getNonLinearCapacityForces(
            self.phys_state.primal_prev,
            self.phys_state.primal_step,
            self.phys_state.time_step,
            self.phys_state.externVar
        )

        EVNL_AS = (
            (1.0 / dt) * self._resi_prev.resi_mass - \
            (1.0 - theta) * self._resi_prev.resi_stress
        )

        residual = Residuals()

        residual.resi_stress = - EVNL_AS + (
            theta * resi_curr.resi_stress +
            (1.0 / dt) * resi_mass
        )

        resi_dual = resi_curr.resi_int - resi_curr.resi_stress

        residual.resi_int = residual.resi_stress + resi_dual

        residual.resi_dual = resi_curr.resi_dual
        residual.resi_ext = theta * resi_curr.resi_ext + (1.0-theta) * self._resi_prev.resi_ext
        residual.resi_cont = resi_curr.resi_cont

        residual.resi = residual.resi_ext - residual.resi_int

        self._temp_stress = stress
        self._temp_internVar = internVar

        resi_curr.resi_mass = resi_mass
        self._resi_temp.append(resi_curr)

        return residual
