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

from ...Objects import DiscreteComputation
from ...Supervis import ConvergenceError
from ...Utilities import PETSc, no_new_attributes, profile
from .iterations_solver import IterationsSolver


def Print(*args):
    print(*args, flush=True)


class SNESSolver(IterationsSolver):
    """Solves a step using PETSc SNES, loops on iterations."""

    _primal_incr = _resi_comp = None
    _scaling = _options = None
    local = snes = fnorm0 = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, local=False):
        super().__init__()
        self._primal_incr = self._resi_comp = None
        self._scaling = self._options = None
        self.local = local
        self.snes = None
        self.fnorm0 = None

    def initialize(self):
        """Initialize the object for the next step."""
        self.current_incr = 0
        # self.current_matrix = None  # keep

    def _evalFunction(self, snes, X, F):
        # Get the solution increment from PETSc
        if snes.getSolutionUpdate().handle:
            self._primal_incr.fromPetsc(snes.getSolutionUpdate(), local=self.local)
        # Increment the solution
        self._primal_incr.applyLagrangeScaling(1 / self._scaling)
        self.state.primal_step += self._primal_incr
        # Build initial residual
        disc_comp = DiscreteComputation(self.problem)
        residual = self.oper.getResidual(self._scaling)
        # Apply Lagrange scaling
        residual.resi.applyLagrangeScaling(1 / self._scaling)
        # Apply DirichletBC into the residual
        diriBCs = disc_comp.getIncrementalDirichletBC(self.state.time_curr, self.state.primal_curr)
        self.current_matrix.applyDirichletBC(diriBCs, residual.resi)
        # Copy to PETSc
        residual.resi.toPetsc(local=self.local).copy(F)
        # print("internVar=", self.oper._tmp_internVar.getValues()[:40:8], flush=True)

    def _evalJacobian(self, snes, X, J, P):
        if self.current_incr % self._matr_update_incr == 0:
            _matrix = self.oper.getJacobian(self.matrix_type)
            self.current_matrix = _matrix
            self._scaling = self.oper.getLagrangeScaling(self.matrix_type)
            _matrix.toPetsc(local=self.local).copy(P)
        self.current_incr += 1
        if J != P:
            J.assemble()  # matrix-free operator
        return PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    def initSNES(self):
        if self.snes and self.local:
            return
        self.oper.initialize()

        self._scaling = self.oper.getLagrangeScaling(self.matrix_type)
        self.current_matrix = self.oper.first_jacobian

        p_jac = self.current_matrix.toPetsc(local=self.local)

        snes = PETSc.SNES().create(comm=p_jac.comm)
        # set a prefix if the solver is local in order to conveniently set options
        if self.local:
            snes.prefix = "lsnes_"
        self.snes = snes

        # register the function in charge of
        # computing the nonlinear residual
        p_resi = p_jac.getVecRight()
        p_resi.set(0)

        self._primal_incr = self.state.primal_step.copy()

        snes.setFunction(self._evalFunction, p_resi)
        snes.setJacobian(self._evalJacobian, p_jac)

        def _monitor(snes, its, fgnorm):
            self.fnorm0 = self.fnorm0 if its else fgnorm
            self.logManager.printConvTableRow(
                [its, fgnorm / self.fnorm0, fgnorm, " - ", self.matrix_type]
            )

        snes.setMonitor(_monitor) if not self.local else None

        rtol = self.get_keyword("CONVERGENCE", "RESI_GLOB_RELA")
        atol = self.get_keyword("CONVERGENCE", "RESI_GLOB_MAXI", 1.0e-24)
        maxiter = self.get_keyword("CONVERGENCE", "ITER_GLOB_MAXI")
        snes.setTolerances(rtol=rtol, atol=atol, max_it=maxiter)

        OptDB = PETSc.Options()
        if not self._options:
            self.linear_solver.build()
            self._options = self.linear_solver.getPetscOptions()
        OptDB.insertString(self._options)
        snes.setFromOptions()

    # @profile
    def solve(self, current_matrix):
        """Solve a step with SNES.

        Raises:
            *ConvergenceError* exception in case of error.
        """

        self.initSNES()

        snes = self.snes
        p_resi, _ = snes.getFunction()

        # solve the nonlinear problem
        b, x = None, p_resi.copy()
        x.set(0)  # zero initial guess
        snes.solve(b, x)

        if snes.getConvergedReason() < 0:
            raise ConvergenceError("MECANONLINE9_7")

        self.oper.finalize()

        return self.current_matrix
