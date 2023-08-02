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

from libaster import deleteTemporaryObjects

from ..Supervis import ConvergenceError
from ..Utilities import with_loglevel, no_new_attributes, profile
from .logging_manager import LoggingManager
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP


class StepSolver(SolverFeature):
    """Solves a step, loops on iterations."""

    provide = SOP.StepSolver
    required_features = [
        SOP.PhysicalProblem,
        SOP.PhysicalState,
        SOP.ConvergenceManager,
        SOP.IncrementalSolver,
    ]

    matr_update_incr = None
    param = None
    current_incr = current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def setParameters(self, param):
        """Set parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param
        assert self._get("NEWTON", "PREDICTION") in (
            "ELASTIQUE",
            "TANGENTE",
        ), f"unsupported value: "

        self.matr_update_incr = self._get("NEWTON", "REAC_ITER", 1)

    def initialize(self):
        """Initialization."""
        self.check_features()
        self.current_incr = -1
        self.current_matrix = None
        self.phys_state.primal_step = self.phys_state.createPrimal(self.phys_pb, 0.0)

        convManager = self.get_feature(SOP.ConvergenceManager)
        convManager.initialize()

        iter_glob = convManager.setdefault("ITER_GLOB_MAXI")
        iter_glob.minValue = 1

    @property
    def contact_manager(self):
        """ContactManager: contact object."""
        return self.get_feature(SOP.Contact, optional=True)

    def update(self, primal_incr, internVar, sigma, convManager):
        """Update the physical state.

        Arguments:
           primal_incr (FieldOnNodes): Displacement increment.
           internVar (FieldOnCells): Internal state variables.
           sigma (FieldOnCells): Stress field.
           convManager (ConvergenceManager): Object that manages the
               convergency criteria.
        """

        self.phys_state.primal_step += primal_incr

        if convManager.isConverged():
            self.phys_state.internVar = internVar
            self.phys_state.stress = sigma

    def _setMatrixType(self):
        """Set matrix type.

        Returns:
            str: Type of matrix to be computed.
        """
        if self.current_incr == 0:
            matrix_type = "PRED_" + self._get("NEWTON", "PREDICTION")
        else:
            matrix_type = self._get("NEWTON", "MATRICE", "TANGENTE")
            if self.current_incr % self.matr_update_incr == 0 or self.contact_manager:
                # make unavailable the current tangent matrix
                self.current_matrix = None
        return matrix_type

    def createLoggingManager(self):
        """Return a logging manager

        Returns:
            LoggingManager: object for logging
        """
        logManager = LoggingManager()
        logManager.addConvTableColumn("NEWTON")
        logManager.addConvTableColumn("RESIDU RELATIF RESI_GLOB_RELA")
        logManager.addConvTableColumn("RESIDU ABSOLU RESI_GLOB_MAXI")
        logManager.addConvTableColumn("RESIDU GEOMETRIQUE RESI_GEOM")
        logManager.addConvTableColumn("OPTION ASSEMBLAGE")

        return logManager

    @profile
    def solve(self, current_matrix):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """

        logManager = self.createLoggingManager()
        logManager.printConvTableEntries()

        self.current_matrix = current_matrix

        convManager = self.get_feature(SOP.ConvergenceManager)
        iter_glob = convManager.setdefault("ITER_GLOB_MAXI")

        incr_solv = self.get_feature(SOP.IncrementalSolver)

        while not convManager.isFinished():
            self.current_incr += 1
            iter_glob.value = self.current_incr

            # Select type of matrix
            matrix_type = self._setMatrixType()

            if self.contact_manager:
                self.contact_manager.pairing(self.phys_pb)

            # Solve current iteration
            primal_incr, internVar, sigma, self.current_matrix = incr_solv.solve(
                matrix_type, self.current_matrix
            )

            # Update
            self.update(primal_incr, internVar, sigma, convManager)

            if self.current_incr > 0:
                logManager.printConvTableRow(
                    [
                        self.current_incr - 1,
                        convManager.get("RESI_GLOB_RELA"),
                        convManager.get("RESI_GLOB_MAXI"),
                        convManager.get("RESI_GEOM"),
                        matrix_type,
                    ]
                )

        if not convManager.isConverged():
            raise ConvergenceError("MECANONLINE9_9")

        deleteTemporaryObjects()

        logManager.printConvTableEnd()

        if self.current_incr % self.matr_update_incr == 0 or self.contact_manager:
            self.current_matrix = None

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        args = self.param
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword])[0].get(parameter, default)

        return args.get(keyword, default)
