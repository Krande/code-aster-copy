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

from ...Supervis import ConvergenceError
from ...Utilities import no_new_attributes, profile
from ..Basics import SolverFeature
from ..Basics import SolverOptions as SOP


class NewtonSolver(SolverFeature):
    """Solves a step, loops on iterations."""

    provide = SOP.ConvergenceCriteria

    required_features = [
        SOP.PhysicalProblem,
        SOP.PhysicalState,
        SOP.ConvergenceManager,
        SOP.IncrementalSolver,
        SOP.LinearSolver,
    ]

    optional_features = [SOP.OperatorsManager, SOP.Contact]

    matr_update_incr = prediction = None
    param = logManager = None
    current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()

    def initialize(self):
        """Initialize the object for the next step."""
        self.check_features()
        self.current_matrix = None
        self.conv_manager.initialize()
        iter_glob = self.conv_manager.setdefault("ITER_GLOB_MAXI")
        iter_glob.minValue = 1

    @property
    def conv_manager(self):
        """ConvergenceManager: Convergence manager object."""
        return self.get_feature(SOP.ConvergenceManager)

    @property
    def contact_manager(self):
        """ContactManager: contact object."""
        return self.get_feature(SOP.Contact, optional=True)

    @property
    def opers_manager(self):
        """OperatorsManager: Operators manager object."""
        return self.get_feature(SOP.OperatorsManager)

    def setParameters(self, param):
        """Assign parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

        self.prediction = self._get("NEWTON", "PREDICTION") or self._get("NEWTON", "MATRICE")

        assert self.prediction in ("ELASTIQUE", "TANGENTE"), f"unsupported value: "

        self.matr_update_incr = self._get("NEWTON", "REAC_ITER", 1)

    def setLoggingManager(self, logManager):
        """Assign the logging manager.

        Arguments:
            logManager (LoggingManager): Logging manager.
        """
        self.logManager = logManager

    def update(self, primal_incr, resi_fields=None, callback=None):
        """Update the physical state.

        Arguments:
            primal_incr (FieldOnNodes): Displacement increment.
            resi_fields (dict of FieldOnNodes): Fields of residual values
        """

        self.phys_state.primal_step += primal_incr

        for key, field in resi_fields.items():
            self.phys_state.set(key, field)

        if callback:
            callback(primal_incr)

    def _resetMatrix(self, current_incr):
        """Reset matrix if needed

        Arguments:
            current_incr (int): index of the current increment

        """
        if (
            self.matr_update_incr > 0 and current_incr % self.matr_update_incr == 0
        ) or self.contact_manager:
            # make unavailable the current tangent matrix
            self.current_matrix = None

    def _setMatrixType(self, current_incr):
        """Set matrix type.

        Arguments:
            current_incr (int): index of the current increment

        Returns:
            str: Type of matrix to be computed.
        """
        if current_incr == 0:
            matrix_type = "PRED_" + self.prediction
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")

            self._resetMatrix(current_incr)
        return matrix_type

    @profile
    def solve(self, current_matrix, callback=None):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        self.current_matrix = current_matrix

        iter_glob = self.conv_manager.setdefault("ITER_GLOB_MAXI")

        incr_solv = self.get_feature(SOP.IncrementalSolver)
        incr_solv.use(self.opers_manager)
        current_incr = -1

        self.opers_manager.initialize()

        while not self.conv_manager.isFinished():
            current_incr += 1

            iter_glob.value = current_incr

            # Select type of matrix
            matrix_type = self._setMatrixType(current_incr)

            # Should the iteration be executed even if the solver converged ?
            force = self.opers_manager.executeIteration(current_incr)

            # Solve current iteration
            primal_incr, self.current_matrix, resi_fields = incr_solv.solve(
                matrix_type, self.current_matrix, force
            )

            # Update
            self.update(primal_incr, resi_fields, callback)

            if current_incr > 0:
                self.logManager.printConvTableRow(
                    [
                        current_incr - 1,
                        self.conv_manager.get("RESI_GLOB_RELA"),
                        self.conv_manager.get("RESI_GLOB_MAXI"),
                        self.conv_manager.get("RESI_GEOM"),
                        matrix_type,
                    ]
                )

        if not self.conv_manager.isConverged():
            raise ConvergenceError("MECANONLINE9_7")

        self.opers_manager.finalize()

        self._resetMatrix(current_incr)

        return self.current_matrix

    def _get(self, keyword, parameter=None, default=None):
        """Return a keyword value"""
        args = self.param
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword][0]).get(parameter, default)

        return args.get(keyword, default)
