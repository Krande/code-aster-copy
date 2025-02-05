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

from abc import ABC, abstractmethod

from ..Basics import SolverOptions as SOP
from ..OperatorsManager import MecaDynaOperatorsManager


class IntegrationType:
    """Time integrators types."""

    Unset = 0x00
    Implicit = 0x01
    Explicit = 0x02
    Multiple = 0x04


class IntegratorName:
    """Time integrators names."""

    Unset = 0x00
    Newmark = 0x01


# FIXME: add ABC after removing SolverFeature
class BaseIntegrator(MecaDynaOperatorsManager):
    """
    Integrator for systems like : M ddX = Fext - Fc(dX) - Fk(X) = funForce(X, dX)
    In case of a linear problem : M ddX = Fext - C dX - K X

    Arguments:
        mass : Mass matrix
        f : difference between external and internal forces
        df : Jacobian matrix of f
    """

    provide = SOP.TimeIntegrator | MecaDynaOperatorsManager.provide

    integration_type = IntegrationType.Unset
    integrator_name = IntegratorName.Unset

    _init_state = _set_up = None

    @classmethod
    def create(cls, name, schema):
        """Setup a solver for the given problem.

        Arguments:
            name : integrator name.
            schema (dict) : *SCHEMA_TEMPS* keyword.

        Returns:
            *BaseIntegrator*: A relevant *BaseIntegrator* object.
        """
        for klass in cls.__subclasses__():
            if klass.integrator_name == name:
                return klass.create(schema)

    def __init__(self):
        super().__init__()
        self._set_up = False
        self._init_state = None

    @property
    def t0(self):
        """Time at the beginning of the step."""
        return self.phys_state.time_prev

    @property
    def dt(self):
        """Returns the current time step"""
        return self.phys_state.time_step

    def getLagrangeScaling(self, matrix_type):
        """Returns Lagrange scaling.

        Arguments:
            matrix_type (str): type of matrix used.

        """

        if self._lagr_scaling is None:
            self._first_jacobian = self.getJacobian(matrix_type)

        return self._lagr_scaling

    def getInitialState(self):
        """Return the physical state used at the beginning of the iteration.

        Returns:
            PhysicalState: Initial physical state.
        """
        return self._init_state

    def setInitialState(self, initial_state):
        """Define the state at the beginning of the iteration.

        Arguments:
            state (PhysicalState): State at the beginning of the iteration.
        """
        self._init_state = initial_state.duplicate()
        if not self._set_up:
            self.setup()

    def initializeStep(self):
        """Define the step parameters."""
        raise NotImplementedError

    def updateVariables(self, q, dq=None, ddq=None):
        raise NotImplementedError

    def getJacobian(self, matrix_type):
        """Compute the jacobian matrix.

        Arguments:
            matrix_type (str): type of matrix used.

        Returns:
            AssemblyMatrixDisplacementReal: Jacobian matrix.
        """
        raise NotImplementedError

    @abstractmethod
    def getResidual(self, scaling=1.0):
        """Compute the residue vector."""
        raise NotImplementedError

    def computeAcceleration(self):
        """Computes the acceleration."""
        force = self.getFunctional(self.t0, self.dt, self.U, self.dU, self.d2U).resi
        linear_solver = self.get_feature(SOP.LinearSolver)
        if not self._mass.isFactorized():
            linear_solver.factorize(self._mass)
        self.phys_state.current.d2U = linear_solver.solve(force)

    @property
    def U0(self):
        """Initial primal unknowns."""
        return self._init_state.current.U

    @U0.setter
    def U0(self, field):
        """Set initial primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self._init_state.current.U = field

    @property
    def dU0(self):
        """Initial derivative of primal unknowns."""
        return self._init_state.current.dU

    @dU0.setter
    def dU0(self, field):
        """Set derivative of initial primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self._init_state.current.dU = field

    @property
    def d2U0(self):
        """Initial second derivative of primal unknowns."""
        return self._init_state.current.d2U

    @d2U0.setter
    def d2U0(self, field):
        """Set second derivative of initial primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self._init_state.current.d2U = field

    @property
    def U(self):
        """Current primal unknowns."""
        return self.phys_state.current.U

    @U.setter
    def U(self, field):
        """Set current primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self.phys_state.current.U = field

    @property
    def dU(self):
        """Current derivative of primal unknowns."""
        return self.phys_state.current.dU

    @dU.setter
    def dU(self, field):
        """Set derivative of current primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self.phys_state.current.dU = field

    @property
    def d2U(self):
        """Current second derivative of primal unknowns."""
        return self.phys_state.current.d2U

    @d2U.setter
    def d2U(self, field):
        """Set second derivative of current primal unknowns.

        Arguments:
            field (FieldOnNodesReal): field value
        """
        self.phys_state.current.d2U = field

    def setup(self):
        """set up the integrator."""
        self._elem_mass = self._getElemMassMatrix()
        self._mass = self._getMassMatrix()
        self._set_up = True
