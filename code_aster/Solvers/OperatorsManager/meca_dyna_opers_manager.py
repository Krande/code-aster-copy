# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
from ...Objects import AssemblyMatrixDisplacementReal, DiscreteComputation
from ...Supervis import IntegrationError

from ..Basics import SolverOptions as SOP


class MecaDynaOperatorsManager(BaseOperatorsManager):
    """Solve an iteration."""

    _mass = _elem_mass = None
    _temp_stress = _temp_internVar = None

    def __init__(self):
        super().__init__()
        self._elem_mass = self._mass = None
        self._temp_stress = self._temp_internVar = None

    def _getElemMassMatrix(self):
        """Compute the elementary mass matrix."""
        disc_comp = DiscreteComputation(self.phys_pb)
        matr_elem_mass = disc_comp.getMassMatrix(
            time=self.phys_state.time_curr, varc_curr=self.phys_state.internVar
        )
        return matr_elem_mass

    def _getMassMatrix(self):
        """Compute the mass matrix."""
        mass_matr = AssemblyMatrixDisplacementReal(self.phys_pb)
        mass_matr.addElementaryMatrix(self._elem_mass)
        mass_matr.assemble()
        return mass_matr

    def getFunctional(self, t, dt, U, dU, d2U, scaling=1.0):
        """Computes the functional."""
        temp_phys_state = self.getTempPhysicalState(t, dt, U, dU, d2U)
        temp_phys_state.swap(self.phys_state)
        resi_state, internVar, stress = super().getResidual(scaling=scaling)
        self.phys_state.swap(temp_phys_state)
        self._temp_stress = stress
        self._temp_internVar = internVar
        return resi_state

    def getStiffAndDamp(self, t, dt, U, dU, d2U, matrix_type):
        """Computes the jacobian."""
        temp_phys_state = self.getTempPhysicalState(t, dt, U, dU, d2U)
        temp_phys_state.swap(self.phys_state)

        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute elementary matrix
        codret, matr_elem_rigi, matr_elem_dual = disc_comp.getInternalTangentMatrix(
            self.phys_state, matrix_type=matrix_type, assemble=False
        )

        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        contact_manager = self.get_feature(SOP.Contact, optional=True)

        matr_elem_cont = disc_comp.getContactTangentMatrix(self.phys_state, contact_manager)

        matr_elem_ext = disc_comp.getExternalTangentMatrix(self.phys_state)

        self.phys_state.swap(temp_phys_state)

        # Assemble stiffness matrix
        rigi_matr = AssemblyMatrixDisplacementReal(self.phys_pb)
        rigi_matr.addElementaryMatrix(matr_elem_rigi)
        rigi_matr.addElementaryMatrix(matr_elem_dual)
        rigi_matr.addElementaryMatrix(matr_elem_cont)
        rigi_matr.addElementaryMatrix(matr_elem_ext)
        rigi_matr.assemble()

        # Assemble damping matrix
        time_curr = t + dt

        matr_elem_damp = disc_comp.getDampingMatrix(
            massMatrix=self._elem_mass,
            stiffnessMatrix=matr_elem_rigi,
            varc_curr=self.phys_state.externVar,
        )

        damp_matr = AssemblyMatrixDisplacementReal(self.phys_pb)
        damp_matr.addElementaryMatrix(matr_elem_damp)
        damp_matr.assemble()

        return rigi_matr, damp_matr

    def getTempPhysicalState(self, t, dt, U, dU, d2U):
        """Creates a temporary physical state"""
        # D'ou vient le primal step qui dois etre utlisée ?
        result = self.phys_state.duplicate()

        result.time_prev = t
        result.time_step = dt
        result.time_curr = t + dt

        result.current.U = U
        result.current.dU = dU
        result.current.d2U = d2U

        return result
