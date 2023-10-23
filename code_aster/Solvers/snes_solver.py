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

from ..Objects import DiscreteComputation
from ..Supervis import ConvergenceError
from ..Utilities import PETSc, no_new_attributes, profile
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP


class SNESSolver(SolverFeature):
    """Solves a step, loops on iterations."""

    provide = SOP.ConvergenceCriteria

    required_features = [SOP.PhysicalProblem, SOP.PhysicalState, SOP.LinearSolver]

    optional_features = [SOP.OperatorsManager]

    param = logManager = None
    matr_update_incr = prediction = None
    current_incr = current_matrix = None
    _primal_incr = _resi_comp = None
    _scaling = _options = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()
        self.matr_update_incr = self.prediction = None
        self.param = self.logManager = None
        self.current_incr = self.current_matrix = None
        self._primal_incr = self._resi_comp = None
        self._scaling = self._options = None

    def initialize(self):
        """Initialize the object for the next step."""
        self.check_features()
        self.current_incr = 0
        self.current_matrix = None

    @property
    def matrix_type(self):
        """Matrix type"""
        return self._setMatrixType()

    @property
    def opers_manager(self):
        """OperatorsManager: Operators manager object."""
        return self.get_feature(SOP.OperatorsManager)

    def setParameters(self, param):
        """Set parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

        self.prediction = self._get("NEWTON", "PREDICTION")

        assert self.prediction in ("ELASTIQUE", "TANGENTE"), f"unsupported value: "

        self.matr_update_incr = self._get("NEWTON", "REAC_ITER", 1)

    def setLoggingManager(self, logManager):
        """Assign the logging manager.

        Arguments:
            logManager (LoggingManager): Logging manager.
        """
        self.logManager = logManager

    def _setMatrixType(self):
        """Set matrix type.

        Returns:
            str: Type of matrix to be computed.
        """
        if self.current_incr == 0:
            matrix_type = "PRED_" + self.prediction
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")

        return matrix_type

    def _evalFunction(self, snes, X, F):
        # Get the solution increment from PETSc
        self._primal_incr.fromPetsc(snes.getSolutionUpdate())
        # Increment the solution
        self._primal_incr.applyLagrangeScaling(1 / self._scaling)
        self.phys_state.primal_step += self._primal_incr
        # Build initial residual
        disc_comp = DiscreteComputation(self.phys_pb)
        residual = self.opers_manager.getResidual(self._scaling)
        # Apply Lagrange scaling
        residual.resi.applyLagrangeScaling(1 / self._scaling)
        # Apply DirichletBC into the residual
        diriBCs = disc_comp.getIncrementalDirichletBC(
            self.phys_state.time_curr, self.phys_state.primal_curr
        )
        self.current_matrix.applyDirichletBC(diriBCs, residual.resi)
        # Copy to PETSc
        residual.resi.toPetsc().copy(F)

    def _evalJacobian(self, snes, X, J, P):
        if self.current_incr % self.matr_update_incr == 0:
            _matrix = self.opers_manager.getJacobian(self.matrix_type)
            self.current_matrix = _matrix
            self._scaling = self.opers_manager.getLagrangeScaling(self.matrix_type)
            _matrix.toPetsc().copy(P)
        self.current_incr += 1
        if J != P:
            J.assemble()  # matrix-free operator
        return PETSc.Mat.Structure.SAME_NONZERO_PATTERN

    @profile
    def solve(self, current_matrix):
        """Solve a step with SNES.

        Raises:
            *ConvergenceError* exception in case of error.
        """

        self.opers_manager.initialize()

        self._scaling = self.opers_manager.getLagrangeScaling(self.matrix_type)
        self.current_matrix = self.opers_manager.first_jacobian

        p_jac = self.current_matrix.toPetsc()

        snes = PETSc.SNES().create()

        # register the function in charge of
        # computing the nonlinear residual
        p_resi = p_jac.getVecRight()
        p_resi.set(0)

        self._primal_incr = self.phys_state.primal_step.copy()

        snes.setFunction(self._evalFunction, p_resi)
        snes.setJacobian(self._evalJacobian, p_jac)

        def _monitor(snes, its, fgnorm):
            self.logManager.printConvTableRow([its, " - ", fgnorm, " - ", self.matrix_type])

        snes.setMonitor(_monitor)

        OptDB = PETSc.Options()
        if not self._options:
            linear_solver = self.get_feature(SOP.LinearSolver)
            linear_solver.build()
            self._options = linear_solver.getPetscOptions()
        OptDB.insertString(self._options)
        snes.setFromOptions()

        # solve the nonlinear problem
        b, x = None, p_resi.copy()
        x.set(0)  # zero initial guess
        snes.solve(b, x)

        if snes.getConvergedReason() < 0:
            raise ConvergenceError("MECANONLINE9_7")

        self.opers_manager.finalize()

        return self.current_matrix

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        args = self.param
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword])[0].get(parameter, default)

        return args.get(keyword, default)
