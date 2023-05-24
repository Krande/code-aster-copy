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

from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP
from ..Supervis import ConvergenceError
from ..Utilities import no_new_attributes, profile


class GeometricSolver(SolverFeature):
    """Solves a step, loops on iterations."""

    provide = SOP.ConvergenceCriteria
    required_features = [
        SOP.PhysicalProblem,
        SOP.PhysicalState,
        SOP.ConvergenceManager,
        SOP.IncrementalSolver,
    ]
    optional_features = [SOP.Contact]

    matr_update_incr = prediction = None
    param = logManager = None
    current_incr = current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()

    def initialize(self):
        """Initialize the object for the next step."""
        self.check_features()
        self.current_incr = 0
        self.current_matrix = None

    @property
    def contact_manager(self):
        """ContactManager: contact object."""
        return self.get_feature(SOP.Contact, optional=True)

    def setParameters(self, param):
        """Assign parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

        self.prediction = self._get("NEWTON", "PREDICTION")
        assert self.prediction in ("ELASTIQUE", "TANGENTE"), f"unsupported value: {self.prediction}"

        self.matr_update_incr = self._get("NEWTON", "REAC_INCR", 1)

    def setLoggingManager(self, logManager):
        """Assign the logging manager.

        Arguments:
            logManager (LoggingManager): Logging manager.
        """
        self.logManager = logManager

    def update(self, primal_incr, internVar, sigma, convManager):
        """Update the physical state.

        Arguments:
            primal_incr (FieldOnNodes): Displacement increment.
            internVar (FieldOnCells): Internal state variables.
            sigma (FieldOnCells): Stress field.
            convManager (ConvergenceManager): Object that manages the
                convergency criteria.
        """
        convManager.evalGeometricResidual(primal_incr)
        self.phys_state.primal_step += primal_incr

        if convManager.hasConverged():
            self.phys_state.internVar = internVar
            self.phys_state.stress = sigma
        elif self._get("CONTACT", "ALGO_RESO_GEOM") == "NEWTON":
            self.contact_manager.update(self.phys_state)
            self.contact_manager.pairing(self.phys_pb)

    def hasFinished(self, convManager):
        """Tell if there are iterations to be computed.

        Arguments:
            convManager (ConvergenceManager): convergence manager.

        Returns:
            bool: *True* if there is no iteration to be computed, *False* otherwise.
        """
        if self.current_incr > self._get("CONVERGENCE", "ITER_GLOB_MAXI"):
            return True
        if self.current_incr < 2:
            return False
        return convManager.hasConverged()

    def _setMatrixType(self):
        """Set matrix type.

        Returns:
            str: Type of matrix to be computed.
        """
        if self.current_incr == 0:
            matrix_type = "PRED_" + self.prediction
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")
            if self.current_incr % self.matr_update_incr == 0 or self.contact_manager:
                # make unavailable the current tangent matrix
                self.current_matrix = None
        return matrix_type

    @profile
    def solve(self, current_matrix):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        self.current_matrix = current_matrix
        convManager = self.get_feature(SOP.ConvergenceManager)
        convManager.initialize()
        iteration = self.get_feature(SOP.IncrementalSolver)

        if self.contact_manager:
            self.contact_manager.pairing(self.phys_pb)

        while not self.hasFinished(convManager):
            # Select type of matrix
            matrix_type = self._setMatrixType()

            # Solve current iteration
            primal_incr, internVar, sigma, self.current_matrix = iteration.solve(
                matrix_type, self.current_matrix
            )

            # Update
            self.update(primal_incr, internVar, sigma, convManager)
            # self.phys_state.debugPrint("<iter+> ")

            if self.current_incr > 0:
                self.logManager.printConvTableRow(
                    [
                        self.current_incr - 1,
                        convManager.getCriteria("RESI_GLOB_RELA"),
                        convManager.getCriteria("RESI_GLOB_MAXI"),
                        convManager.getCriteria("RESI_GEOM"),
                        matrix_pred,
                    ]
                )

            self.current_incr += 1
            matrix_pred = matrix_type

        if not convManager.hasConverged():
            raise ConvergenceError("MECANONLINE9_7")

        # print(f"| Nombre d'itérations de Newton : {self.current_incr - 1}")
        if self.current_incr % self.matr_update_incr == 0 or self.contact_manager:
            self.current_matrix = None

        return self.current_matrix

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        args = self.param
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword])[0].get(parameter, default)

        return args.get(keyword, default)
