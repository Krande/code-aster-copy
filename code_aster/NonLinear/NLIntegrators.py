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

#!/usr/bin/env python3

import unittest
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
from code_aster.NonLinear.nl_integr_features import NLIntegrFeature
from code_aster.NonLinear.nl_integr_features import NLIntegrOptions as FOP
from time_stepper import TimeStepper

PLOT = False

# already exists in base_utils
def force_list(values):
    """Ensure `values` is iterable (list, tuple, array...) and return it as
    a list."""
    if not isinstance(values, (list, tuple)):
        values = [values]
    return list(values)


class PhysicalState:
    """This is the physical state at one time step.

    Arguments:
        U (*FieldOnNodesReal*, optional): Initial field of primal unknowns.
        dU (*FieldOnNodesReal*, optional): Initial field of derivative of
            primal unknowns.
        d2U (*FieldOnNodesReal*, optional): Initial field of second
            derivative of primal unknowns.
    """

    def __init__(self, U=None, dU=None, d2U=None):
        super().__init__()
        # self._time = time ?
        self._primal = U
        self._dprim = dU
        self._d2prim = d2U

    def copy(self):
        """Return a copy.

        Returns:
            PhysicalState: New object.
        """
        return PhysicalState(self.U.copy(), self.dU.copy(), self.d2U.copy())

    @property
    def U(self):
        """FieldOnNodesReal: Primal field."""
        return self._primal

    @U.setter
    def U(self, field):
        self._primal = field

    @property
    def dU(self):
        """FieldOnNodesReal: Primal field."""
        return self._dprim

    @dU.setter
    def dU(self, field):
        self._dprim = field

    @property
    def d2U(self):
        """FieldOnNodesReal: Primal field."""
        return self._d2prim

    @d2U.setter
    def d2U(self, field):
        self._d2prim = field

    # FIXME tmp: for unittests
    def tolist(self):
        return self._primal, self._dprim, self._d2prim

    def __repr__(self):
        """Representation of the values."""
        return f"({self.U}, {self.dU}, {self.d2U})"


class StorageManager(NLIntegrFeature):
    """Object that stores the states of the problem."""

    provide = FOP.Storage

    def __init__(self):
        super().__init__()
        self._result = {}

    def add(self, time, state):
        """Store the current state.

        Arguments:
            time (float): Current state.
            state (PhysicalState): Current physical state.
        """
        assert time not in self._result, f"already exists: {time:.3f}, values: {self._result[time]}"
        assert isinstance(state, PhysicalState), f"{type(state)} {state!r}"
        self._result[time] = state

    def get(self, time=None):
        """Return the physical state at the given time.

        Arguments:
            time (float, optional): Time of extraction. Use the last time if
                it is not provided.

        Returns:
            PhysicalState: Physical state.
        """
        return self._result[time if time is not None else self.getLastTime()]

    def getLastTime(self):
        """Return the last time already stored.

        Returns:
            float: Last time.
        """
        assert self._result, "expecting at least the initial state"
        return max(self._result.keys())

    def initialize(self, time, initial_state):
        """Define the initial state.

        Arguments:
            time (float): Initial time.
            initial_state (*PhysicalState*): initial state.
        """
        self.add(time, initial_state)


class TransientHistory(NLIntegrFeature):
    """Solve a transient problem over a list of times."""

    provide = FOP.ProblemSolver
    required_features = [FOP.TimeStepper, FOP.Storage, FOP.StepSolver]

    def __init__(self):
        super().__init__()

    def initialize(self, initial_time, initial_state):
        """Define the transient.

        Arguments:
            initial_time (float): Initial time.
            initial_state (PhysicalState): State at the initial time.
        """
        self.check_features()
        steps = self.get_feature(FOP.TimeStepper)
        storage = self.get_feature(FOP.Storage)
        steps.setInitialStep(initial_time)
        storage.initialize(initial_time, initial_state)

    def solve(self, start=None, end=None):
        """Solve the system from *starttime* to *end*.

        Arguments:
            start (float, optional): Start solving at this time.
            end (float, optional): Last time to be solved.
        """
        self.check_features()
        storage = self.get_feature(FOP.Storage)
        last = storage.getLastTime()
        if start is not None:
            last = max(last, start)
        steps = self.get_feature(FOP.TimeStepper)
        steps.setInitialStep(last)
        steps.setFinalStep(end)

        solver = self.get_feature(FOP.StepSolver)
        # Main loop
        current = steps.getCurrent()
        steps.start()
        while not steps.hasFinished():
            solver.setInitialState(storage.get(current))
            next_time = steps.getCurrent()
            delta_t = steps.getIncrement()
            print("*" * 70)
            print("Solving from {} to {}".format(current, current + delta_t))
            solver.solve(current, delta_t)
            storage.add(next_time, solver.getCurrentState())
            steps.completed()
            current = next_time

    # FIXME tmp: for the unittests
    @property
    def solutions(self):
        storage = self.get_feature(FOP.Storage)
        return storage._result


class Options:
    """Integration options."""

    Unset = 0x00
    Implicit = 0x01
    Explicit = 0x02
    Multiple = 0x04


class StepSolver(NLIntegrFeature):
    """Solve a step with given integrator, loops on iterations."""

    provide = FOP.StepSolver
    options = Options.Unset

    @classmethod
    def setup(cls, integrator, coef=None):
        """Setup a solver for the given integrator.

        Arguments:
            integrator (list[Integrator]): one or more *Integrator* objects.
            coef (list[float]): Splitting coefficients if more than one
                integrator (size = number of integrators - 1).

        Returns:
            *StepSolver*: A relevant *StepSolver* object.
        """

        def _setup_one(integrator):
            found = None
            for klass in cls.__subclasses__():
                if klass.support(integrator):
                    found = klass(integrator)
                    break
            assert found, f"unsupported integrator: {integrator}"
            return found

        integrator = force_list(integrator)
        multi = len(integrator) > 1
        if multi:
            assert coef is not None
            coef = force_list(coef)
            assert len(integrator) == len(coef) + 1
            assert len(coef) == 1, "only two steps supported yet"
            solver = MultiStepsSolver([_setup_one(i) for i in integrator], coef)
        else:
            solver = _setup_one(integrator[0])
        return solver

    def __init__(self, integrator, rtol=1.0e-6, atol=1.0e-12, maxits=200):
        super().__init__()
        self._resi0 = None
        self._integrator = integrator
        self._rtol = rtol
        self._atol = atol
        self._maxits = maxits

    def setInitialState(self, initial_state):
        """Define the initial state of the integrator.

        Arguments:
            initial_state (PhysicalState): State at the beginning of the iteration.
        """
        return self._integrator.setInitialState(initial_state)

    def setIntermediateState(self, state):
        """Define the state at the beginning of the iteration.

        Arguments:
            state (PhysicalState): State at the beginning of the sub-step.
        """
        assert hasattr(self._integrator, "setIntermediateState"), "only for OnSubStepIntegrator"
        return self._integrator.setIntermediateState(state)

    def getInitialState(self):
        """Return the physical state used at the beginning of the iteration.

        Returns:
            PhysicalState: Initial physical state.
        """
        return self._integrator.getInitialState()

    def getCurrentState(self):
        """Return the state at the end of the step.

        Returns:
            tuple[PhysicalState]: Current state (as tuple).
        """
        return self._integrator.getCurrentState()

    def solve(self, t_init, delta_t):
        """Solve a step.

        Arguments:
            t_init (float): Time at the beginning of the step.
            delta_t (float): Time step.
        """
        raise NotImplementedError()

    @classmethod
    def support(cls, integrator):
        """Tell if the StepSolver supports this kind of integrator.

        Arguments:
            integrator (Integrator): *Integrator* object.

        Returns:
            bool: *True* if the integrator is supported, *False* otherwise.
        """
        return cls.options == integrator.options


class MultiStepsSolver(StepSolver):
    """Solve a step using two solvers and integrators.

    Arguments :
        stepSolvers (list[StepSolver]): step solvers
        coef (list[float]): list of coefficients to subdivise the time step:
            [t_init, t_init + c0 * delta_t, t_init + delta_t]
    """

    options = Options.Multiple

    def __init__(self, step_solvers: list, coef: list):
        self._resi0 = None
        assert len(step_solvers) == 2, len(step_solvers)
        assert isinstance(step_solvers[0]._integrator, Integrator)
        assert isinstance(step_solvers[1]._integrator, OnSubStepIntegrator)
        self._step_solvers = step_solvers
        assert len(coef) == 1, "only two substeps supported yet"
        self._coef = coef[0]

    def setInitialState(self, initial_state):
        """Define the initial state of the integrator.

        Arguments:
            state (PhysicalState): State at the beginning of the iteration.
        """
        self._step_solvers[0].setInitialState(initial_state)

    def getCurrentState(self):
        """Return the state at the end of the step."""
        return self._step_solvers[-1].getCurrentState()

    def solve(self, t_init, delta_t):
        step0, step1 = self._step_solvers
        print("++ Solving stage 1")
        state0 = step0.getInitialState()
        step0.solve(t_init, self._coef * delta_t)
        state_interm = step0.getCurrentState().copy()
        print("++ Solving stage 2")
        step1.setInitialState(state0)
        step1.setIntermediateState(state_interm)
        step1.solve(t_init, delta_t)


class ImplicitStepSolver(StepSolver):
    """Solve a step, loop on iterations."""

    options = Options.Implicit

    def hasConverged(self, resi, iterat):
        """Tell if the step has converged.

        Arguments:
            resi (float): Computed residue.
            iterat (int): Current iteration index.

        Returns:
            bool: Convergence state.
        """
        if iterat > self._maxits:
            raise RuntimeError(f"ConvergenceError: Could not integrate in {iterat} iterations")
        return (
            np.linalg.norm(resi) < self._rtol * np.linalg.norm(self._resi0)
            or np.linalg.norm(resi) < self._atol
        )

    def solve(self, t_init, delta_t):
        """Solve a step.

        Arguments:
            t_init (float): Time at the beginning of the step.
            delta_t (float): Time step.
        """
        iterat = -1
        intg = self._integrator
        intg.initializeStep(t_init, delta_t)
        jac = intg.computeJacobian()
        self._resi0 = intg.computeResidue()
        resi = self._resi0
        while not self.hasConverged(resi, iterat + 1):
            iterat += 1
            if iterat > 0:
                jac = intg.computeJacobian()
                resi = intg.computeResidue()
            print(f"-- iter = {iterat:3d}, ||err|| = {np.linalg.norm(resi):g}")
            deltaU = np.linalg.solve(jac, -resi)
            intg.updateVariables(deltaU)

        intg.setInitialState(intg.getCurrentState())


class ExplicitStepSolver(StepSolver):
    """Solves a step, loops on iterations."""

    options = Options.Explicit

    def solve(self, t_init, delta_t, *args):
        """Solve a step.

        Arguments:
            t_init (float): Time at the beginning of the step.
            delta_t (float): Time step.
            args (misc): Additional arguments passed to the integrator.
        """
        self._integrator.initializeStep(t_init, delta_t)
        self._integrator.integrate()


class Integrator:
    """
    Integrator for systems like : M ddX = Fext - Fc(dX) - Fk(X) = funForce(X, dX)
    In case of a linear problem : M ddX = Fext - C dX - K X
    Arguments :
        mass : Mass matrix
        f : difference between external and internal forces
        df : Jacobian matrix of f
    """

    options = Options.Unset

    def __init__(self, mass, f, df):
        self._t0 = None
        self._dt = None
        self._mass = mass
        self._funcF = f
        self._funcDF = df
        self._state0 = None
        self._state = None

    def getInitialState(self):
        """Return the physical state used at the beginning of the iteration.

        Returns:
            PhysicalState: Initial physical state.
        """
        return self._state0

    def getCurrentState(self):
        """Return the physical state at the current iteration.

        Returns:
            PhysicalState: Current physical state (at the end of the iteration).
        """
        return self._state

    def setInitialState(self, initial_state):
        """Define the state at the beginning of the iteration.

        Arguments:
            state (PhysicalState): State at the beginning of the iteration.
        """
        self._state0 = initial_state.copy()

    def initializeStep(self, t_init, delta_t):
        """Define the step parameters.

        Arguments:
            t_init (float): Time at the beginning of the step.
            delta_t (float): Time step.
        """
        self._t0 = t_init
        self._dt = delta_t

    def updateVariables(self, q, dq=None, ddq=None):
        raise NotImplementedError()

    def computeJacobian(self):
        """Compute the jacobian matrix."""
        raise NotImplementedError()

    def computeResidue(self):
        """Compute the residue vector."""
        raise NotImplementedError()

    def computeAcceleration(self):
        force = self._funcF(self._t0, self._dt, self.U, self.dU, self.d2U)
        self._state.d2U = np.linalg.solve(self._mass, force)

    @property
    def U(self):
        """Current primal unknowns."""
        return self._state.U

    @property
    def U0(self):
        """Initial primal unknowns."""
        return self._state0.U

    @property
    def dU(self):
        """Current derivative of primal unknowns."""
        return self._state.dU

    @property
    def dU0(self):
        """Initial derivative of primal unknowns."""
        return self._state0.dU

    @property
    def d2U(self):
        """Current second derivative of primal unknowns."""
        return self._state.d2U

    @property
    def d2U0(self):
        """Initial second derivative of primal unknowns."""
        return self._state0.d2U

    def __repr__(self):
        """Representation of an Integrator with its states values."""
        string = f"{self.__class__.__name__}"  # 0x{hex(id(self))}
        if self._state0:
            string += f" init: {self._state0!r}"
        if self._state:
            string += f" curr: {self._state!r}"
        if hasattr(self, "_state_m") and self._state_m:
            string += f" second: {self._state_m!r}"
        return string


class OnSubStepIntegrator(Integrator):
    """
    Integrator that works on a sub-step.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._state_m = None

    def setIntermediateState(self, state):
        """Define the state at the beginning of the iteration.

        Arguments:
            state (PhysicalState): State at the beginning of the sub-step.
        """
        self._state_m = state.copy()


class TRIntegrator(Integrator):
    """
    Implementation of Trapezoidal rule scheme
    """

    options = Options.Implicit

    def __init__(self, mass, f, df):
        super().__init__(mass, f, df)
        self._G0 = np.zeros((mass.shape[0],))
        self._dFc = None

    def initializeStep(self, t_init, delta_t):
        """Define the step parameters.

        Arguments:
            t_init (float): Time at the beginning of the step.
            delta_t (float): Time step.
        """
        super().initializeStep(t_init, delta_t)
        self._state = self._state0.copy()
        self.computeAcceleration()
        # print("++ init", repr(self))

    def computeJacobian(self):
        """Compute the jacobian matrix."""
        dt = self._dt
        dFk, dFc = self._funcDF(self._t0, self._dt, self.U, self.dU, self.d2U)
        K = -dFk
        C = -dFc
        self._dFc = dFc
        jac = self._mass * 4.0 / dt**2 + 2.0 * C / dt + K
        return jac

    def computeResidue(self):
        """Compute the residue vector."""
        dt = self._dt
        G0 = self.U - self.U0 - 0.5 * dt * (self.dU + self.dU0)
        G1 = self.dU - self.dU0 - 0.5 * dt * (self.d2U + self.d2U0)

        residue = self._mass.dot(G1) * 2 / dt + (self._mass * 4 / dt**2 + 2 / dt * self._dFc).dot(
            G0
        )
        self._G0 = G0
        return residue

    def updateVariables(self, deltaU):
        """Update the physical state."""
        dt = self._dt
        self._state.U += deltaU
        self._state.dU += 2 * (deltaU + self._G0) / dt
        self.computeAcceleration()


class NewmarkIntegrator(Integrator):
    """
    Implementation of Newmark's scheme
    """

    options = Options.Implicit

    def __init__(self, mass, funForce, funDF, gamma=0.5, beta=0.25):
        super().__init__(mass, funForce, funDF)
        self._gamma = gamma
        self._beta = beta
        self._J = np.zeros_like(mass)

    def initializeStep(self, t_init, delta_t):
        """Define the step parameters.

        Arguments:
            t_init (float): Time at the beginning of the step.
            delta_t (float): Time step.
        """
        super().initializeStep(t_init, delta_t)
        dt = self._dt
        self._state = PhysicalState()
        self._state.U = self.U0 + dt * self.dU0 + (0.5 - self._beta) * dt**2 * self.d2U0
        self._state.dU = self.dU0 + (1 - self._gamma) * dt * self.d2U0
        self._state.d2U = np.zeros_like(self.U)

    def computeJacobian(self):
        """Compute the jacobian matrix."""
        dt = self._dt
        beta = self._beta
        gamma = self._gamma
        dFk, dFc = self._funcDF(self._t0, self._dt, self.U, self.dU, self.d2U)
        K = -dFk
        C = -dFc
        self._dFc = dFc
        jac = K + gamma / beta / dt * C + 1 / beta / dt**2 * self._mass
        return jac

    def computeResidue(self):
        """Compute the residue vector."""
        force = self._funcF(self._t0, self._dt, self.U, self.dU, self.d2U)
        residue = self._mass.dot(self.d2U) - force
        return residue

    def updateVariables(self, deltaU):
        """Update the physical state."""
        dt = self._dt
        self._state.U += deltaU
        self._state.dU += self._gamma / self._beta / dt * deltaU
        self._state.d2U += 1 / self._beta / dt**2 * deltaU


class BDF2Integrator(OnSubStepIntegrator):
    """
    Implementation of Backward Differentiation Formula scheme of order 2
    """

    options = Options.Implicit

    def __init__(self, mass, funForce, funDF, gamma=2 - sqrt(2)):
        super().__init__(mass, funForce, funDF)
        self._gamma = gamma
        self._G0 = np.zeros_like(mass)

    def initializeStep(self, t_init, delta_t):
        """Define the step parameters.

        Arguments:
            t_init (float): Time at the beginning of the step.
            delta_t (float): Time step.
        """
        super().initializeStep(t_init, delta_t)
        self._state = self._state0.copy()
        self.computeAcceleration()
        # print("++ init", repr(self))

    def computeJacobian(self):
        """Compute the jacobian matrix."""
        dt = self._dt
        g = self._gamma
        dFk, dFc = self._funcDF(self._t0, self._dt, self.U, self.dU, self.d2U)
        K = -dFk
        C = -dFc
        self._dFc = dFc
        g3 = (1 - g) / (2 - g)
        jac = self._mass / g3**2 / dt**2 + C / g3 / dt + K
        return jac

    def computeResidue(self):
        """Compute the residue vector."""
        dt = self._dt
        g = self._gamma
        g1 = 1 / g / (2 - g)
        g2 = -((1 - g) ** 2) / g / (2 - g)
        g3 = (1 - g) / (2 - g)
        G0 = self.U - g1 * self.Um - g2 * self.U0 - g3 * dt * self.dU
        self._G0 = G0
        G1 = self.dU - g1 * self.dUm - g2 * self.dU0 - g3 * dt * self.d2U
        residue = self._mass.dot(G1) / g3 / dt + (
            self._mass / g3**2 / dt**2 + self._dFc / g3 / dt
        ).dot(G0)
        return residue

    def updateVariables(self, deltaU):
        """Update the physical state."""
        dt = self._dt
        g = self._gamma
        g3 = (1 - g) / (2 - g)
        self._state.U += deltaU
        self._state.dU += (deltaU + self._G0) / g3 / dt
        self.computeAcceleration()

    @property
    def Um(self):
        """Initial primal unknowns."""
        return self._state_m.U

    @property
    def dUm(self):
        """Initial derivative of primal unknowns."""
        return self._state_m.dU

    @property
    def d2Um(self):
        """Initial second derivative of primal unknowns."""
        return self._state_m.d2U


class RK4Integrator(Integrator):
    """
    Implementation of Runge-Kutta order 4 scheme
    """

    options = Options.Explicit

    def __init__(self, masse, funForce, funDF):
        super().__init__(masse, funForce, funDF)
        dim = self._mass.shape[0]
        self._k1 = np.zeros((dim,))
        self._k2 = np.zeros((dim,))
        self._k3 = np.zeros((dim,))
        self._k4 = np.zeros((dim,))
        self._first = True

    def integrate(self):
        """Explicit integration."""
        if self._first:
            self._first = False
            self.dU0[:] = 0.0
            self.d2U0[:] = 0.0

        self._state = self._state0.copy()

        t_init = self._t0
        dt = self._dt

        self._k1 = self.d2U0
        self._k2 = self._funcF(
            t_init,
            0.5 * dt,
            self.U0 + 0.5 * dt * self.dU0,
            self.dU0 + 0.5 * dt * self._k1,
            self.d2U0,
        )
        self._k3 = self._funcF(
            t_init,
            0.5 * dt,
            self.U0 + 0.5 * dt * self.dU0 + 0.25 * dt**2 * self._k1,
            self.dU0 + 0.5 * dt * self._k2,
            self.d2U0,
        )
        self._k4 = self._funcF(
            t_init,
            dt,
            self.U0 + dt * self.dU0 + 0.5 * dt**2 * self._k2,
            self.dU0 + dt * self._k3,
            self.d2U0,
        )

        self._state.U = self.U0 + dt * self.dU0 + dt**2 * (self._k1 + self._k2 + self._k3) / 6.0
        self._state.dU = (
            self.dU0 + dt * (self._k1 + 2.0 * self._k2 + 2.0 * self._k3 + self._k4) / 6.0
        )
        self._state.d2U = self._funcF(t_init, dt, self.U, self.dU, self.d2U0)
        self._state0 = self._state.copy()
        return True


# --- test that the system works


class VanDerPolTest(unittest.TestCase):
    def _buildSystem(self):
        M = np.diag([1.0])

        def funForce(t_init, delta_t, q, dq, ddq=None):
            return -q - (q**2 - 6) * dq

        def df(t_init, delta_t, q, dq, ddq=None):
            return -np.diag([1.0]), np.zeros((1, 1))

        tmax = 1
        ref_tmax = (
            np.array([-3.8851826503372524]),
            np.array([-9.355367190054261]),
            np.array([88.96886344061366]),
        )

        return M, funForce, df, tmax, ref_tmax

    def plotSolution(self, timespan, sol, ref=None):

        plt.rcParams.update({"font.size": 14})
        plt.ioff()
        plt.figure()
        plt.xlabel("time")
        plt.ylabel("signal")
        for dof in range(len(sol[0].tolist()[0])):
            comp = 0
            plt.plot(timespan, [sol[t].tolist()[comp][dof] for t in timespan], "-o", label="$X$ ")
            if ref:
                plt.plot(timespan, ref(timespan)[comp][dof], "-", label="$X$ Ref")
            plt.grid()
            plt.legend()
            plt.show()

            plt.figure()
            plt.xlabel("time")
            plt.ylabel("signal")
            comp = 1
            plt.plot(
                timespan, [sol[t].tolist()[comp][dof] for t in timespan], "-o", label="$\dot{X}$"
            )
            if ref:
                plt.plot(timespan, ref(timespan)[comp][dof], "-", label="$\dot{X}$ Ref")
            plt.grid()
            plt.legend()
            plt.show()

            plt.figure()
            plt.xlabel("time")
            plt.ylabel("signal")
            comp = 2
            plt.plot(
                timespan, [sol[t].tolist()[comp][dof] for t in timespan], "-o", label="$\ddot{X}$"
            )
            if ref:
                plt.plot(timespan, ref(timespan)[comp][dof], "-", label="$\ddot{X}$ Ref")
            plt.grid()
            plt.legend()
            plt.show()

    def testTR(self):
        M, funForce, df, tmax, ref_tmax = self._buildSystem()
        npas = 40
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        integrator = TRIntegrator(M, funForce, df)
        stepsolver = StepSolver.setup(integrator)

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(0.0, PhysicalState(np.array([1.0]), np.array([0.0]), np.array([0.0])))
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions)

        _ = [
            self.assertAlmostEqual(v1, v2, delta=2.0e-2)
            for t1, t2 in zip(ref_tmax, ts.solutions[tmax].tolist())
            for v1, v2 in zip(t1, t2)
        ]

    def testNewmark(self):
        M, funForce, df, tmax, ref_tmax = self._buildSystem()
        npas = 40
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        integrator = NewmarkIntegrator(M, funForce, df)
        stepsolver = StepSolver.setup(integrator)

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(0.0, PhysicalState(np.array([1.0]), np.array([0.0]), np.array([0.0])))
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions)

        _ = [
            self.assertAlmostEqual(v1, v2, delta=2.0e-2)
            for t1, t2 in zip(ref_tmax, ts.solutions[tmax].tolist())
            for v1, v2 in zip(t1, t2)
        ]

    def testRK4(self):
        M, funForce, df, tmax, ref_tmax = self._buildSystem()
        npas = 4000
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        integrator = RK4Integrator(M, funForce, df)
        stepsolver = StepSolver.setup(integrator)

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(0.0, PhysicalState(np.array([1.0]), np.array([0.0]), np.array([0.0])))
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions)

        _ = [
            self.assertAlmostEqual(v1, v2, delta=3)
            for t1, t2 in zip(ref_tmax, ts.solutions[tmax].tolist())
            for v1, v2 in zip(t1, t2)
        ]

    def testTRBDF2(self):
        M, funForce, df, tmax, ref_tmax = self._buildSystem()
        npas = 40
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        stepsolver = StepSolver.setup(
            (TRIntegrator(M, funForce, df), BDF2Integrator(M, funForce, df)), coef=2.0 - sqrt(2)
        )
        # stepsolver = StepSolver()
        # stepsolver.use(TRIntegrator(M, funForce, df))
        # stepsolver.use(BDF2Integrator(M, funForce, df))
        # stepsolver.use(Splitting(2.0 - sqrt(2)))

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(0.0, PhysicalState(np.array([1.0]), np.array([0.0]), np.array([0.0])))
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions)

        _ = [
            self.assertAlmostEqual(v1, v2, delta=4.5e-1)
            for t1, t2 in zip(ref_tmax, ts.solutions[tmax].tolist())
            for v1, v2 in zip(t1, t2)
        ]


class TwoMassesStiffTest(unittest.TestCase):
    def _buildSystem(self):
        M = np.diag([1.0, 1.0])
        k1, k2 = 1e7, 1.0
        K = np.zeros((2, 2))
        K[0, 0] = k1 + k2
        K[0, 1] = -k2
        K[1, 0] = -k2
        K[1, 1] = k2
        wp = 2

        def funForce(t_init, delta_t, q, dq, ddq=None):
            """f = Fext - K q"""
            resu = -np.dot(K, q)
            resu[0] += k1 * np.sin(wp * (t_init + delta_t))
            return resu

        def df(t_init, delta_t, q, dq, ddq=None):
            return -1.0 * K, np.zeros((2, 2))

        tmax = 5

        # Analytical solution
        def ref(t):
            return (
                np.array([np.sin(wp * t), 1.0 / 3 * (wp * np.sin(t) - np.sin(wp * t))]),
                np.array([wp * np.cos(wp * t), 1.0 / 3 * (wp * np.cos(t) - wp * np.cos(2 * t))]),
                np.array(
                    [
                        -(wp**2) * np.sin(wp * t),
                        1.0 / 3 * (-wp * np.sin(t) + wp**2 * np.sin(2 * t)),
                    ]
                ),
            )

        return M, funForce, df, tmax, ref

    def plotSolution(self, timespan, sol, ref=None):

        plt.rcParams.update({"font.size": 14})
        plt.ioff()
        plt.figure()
        plt.xlabel("time")
        plt.ylabel("signal")
        for dof in range(len(sol[0].tolist()[0])):
            comp = 0
            plt.plot(timespan, [sol[t].tolist()[comp][dof] for t in timespan], "-o", label="$X$ ")
            if ref:
                plt.plot(timespan, ref(timespan)[comp][dof], "-", label="$X$ Ref")
            plt.grid()
            plt.legend()
            plt.show()

            plt.figure()
            plt.xlabel("time")
            plt.ylabel("signal")
            comp = 1
            plt.plot(
                timespan, [sol[t].tolist()[comp][dof] for t in timespan], "-o", label="$\dot{X}$"
            )
            if ref:
                plt.plot(timespan, ref(timespan)[comp][dof], "-", label="$\dot{X}$ Ref")
            plt.grid()
            plt.legend()
            plt.show()

            plt.figure()
            plt.xlabel("time")
            plt.ylabel("signal")
            comp = 2
            plt.plot(
                timespan, [sol[t].tolist()[comp][dof] for t in timespan], "-o", label="$\ddot{X}$"
            )
            if ref:
                plt.plot(timespan, ref(timespan)[comp][dof], "-", label="$\ddot{X}$ Ref")
            plt.grid()
            plt.legend()
            plt.show()

    def testTR(self):
        M, funForce, df, tmax, ref = self._buildSystem()
        npas = 40
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        integrator = TRIntegrator(M, funForce, df)
        stepsolver = StepSolver.setup(integrator)

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(
            0.0, PhysicalState(np.array([0.0, 0.0]), np.array([2.0, 0.0]), np.array([0.0, 0.0]))
        )
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions, ref=ref)

        # the acceleration is not tested
        _ = [
            self.assertAlmostEqual(v1, v2, delta=2.0e-2)
            for t1, t2 in zip(ref(tmax)[0:1], ts.solutions[tmax].tolist()[0:1])
            for v1, v2 in zip(t1, t2)
        ]

    def testNewmark(self):
        M, funForce, df, tmax, ref = self._buildSystem()
        npas = 40
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        integrator = NewmarkIntegrator(M, funForce, df)
        stepsolver = StepSolver.setup(integrator)

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(
            0.0, PhysicalState(np.array([0.0, 0.0]), np.array([2.0, 0.0]), np.array([0.0, 0.0]))
        )
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions, ref=ref)

        # the acceleration is not tested
        _ = [
            self.assertAlmostEqual(v1, v2, delta=2.0e-2)
            for t1, t2 in zip(ref(tmax)[0:1], ts.solutions[tmax].tolist()[0:1])
            for v1, v2 in zip(t1, t2)
        ]

    def testRK4(self):
        M, funForce, df, tmax, ref = self._buildSystem()
        npas = 7000
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        integrator = RK4Integrator(M, funForce, df)
        stepsolver = StepSolver.setup(integrator)

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(
            0.0, PhysicalState(np.array([0.0, 0.0]), np.array([2.0, 0.0]), np.array([0.0, 0.0]))
        )
        ts.solve(end=timespan[npas // 3])
        ts.solve(end=timespan[npas // 3] * 2)
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions, ref=ref)

        _ = [
            self.assertAlmostEqual(v1, v2, delta=2.0e-1)
            for t1, t2 in zip(ref(tmax), ts.solutions[tmax].tolist())
            for v1, v2 in zip(t1, t2)
        ]

    def testTRBDF2(self):
        M, funForce, df, tmax, ref = self._buildSystem()

        npas = 80
        timespan = np.linspace(0, tmax, npas + 1)

        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        stepsolver = StepSolver.setup(
            (TRIntegrator(M, funForce, df), BDF2Integrator(M, funForce, df)), coef=2.0 - sqrt(2)
        )

        ts = TransientHistory()
        ts.use(timestepper)
        ts.use(storage)
        ts.use(stepsolver)
        ts.initialize(
            0.0, PhysicalState(np.array([0.0, 0.0]), np.array([2.0, 0.0]), np.array([0.0, 0.0]))
        )
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions, ref=ref)

        _ = [
            self.assertAlmostEqual(v1, v2, delta=2.0e-1)
            for t1, t2 in zip(ref(tmax), ts.solutions[tmax].tolist())
            for v1, v2 in zip(t1, t2)
        ]

    def testTR_TRBDF2(self):
        M, funForce, df, tmax, ref = self._buildSystem()

        npas = 80
        timespan = np.linspace(0, tmax, npas + 1)

        integ1 = TRIntegrator(M, funForce, df)
        integ2 = BDF2Integrator(M, funForce, df)
        timestepper = TimeStepper(timespan)
        storage = StorageManager()
        stepsolver1 = StepSolver.setup(integ1)
        stepsolver2 = StepSolver.setup((integ1, integ2), coef=2.0 - sqrt(2))

        # WARNING: using Trapezoidal scheme for the first 2 steps is not
        # as accurated as with TR+BDF2 but sufficient to check the requirement.
        ts = TransientHistory()
        ts.use(timestepper)
        # we could start with one integrator...
        ts.use(storage)
        ts.use(stepsolver1)
        ts.initialize(
            0.0, PhysicalState(np.array([0.0, 0.0]), np.array([2.0, 0.0]), np.array([0.0, 0.0]))
        )
        ts.solve(0.0, 2.0 * tmax / npas + 1.0e-6)
        # ... and continue with anothers
        ts.discard(stepsolver1)
        ts.use(stepsolver2)
        ts.solve()

        if PLOT:
            self.plotSolution(timespan, ts.solutions, ref=ref)

        _ = [
            self.assertAlmostEqual(v1, v2, delta=2.0e-1)
            for t1, t2 in zip(ref(tmax), ts.solutions[tmax].tolist())
            for v1, v2 in zip(t1, t2)
        ]


class SystemTest(unittest.TestCase):
    """Check for internal methods."""

    def test_options(self):
        M = np.diag([1.0, 1.0])

        integr = TRIntegrator(M, None, None)
        self.assertEqual(integr.options, Options.Implicit)
        solver = StepSolver.setup(TRIntegrator)
        self.assertIsInstance(solver, ImplicitStepSolver)

        solver = StepSolver.setup(RK4Integrator)
        self.assertIsInstance(solver, ExplicitStepSolver)

        class BadIntegrator(Integrator):
            pass

        with self.assertRaisesRegex(AssertionError, "unsupported"):
            StepSolver.setup(BadIntegrator)


if __name__ == "__main__":
    unittest.main()
