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
        SOP.ConvergenceCriteria,
    ]

    current_incr = param = None
    geom = geom_step = current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()

    def setParameters(self, param):
        """Set parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

    def initialize(self):
        """Initialization."""
        self.check_features()
        self.current_incr = -1
        self.phys_state.primal_step = self.phys_state.createPrimal(self.phys_pb, 0.0)
        self.geom_step = self.phys_state.createPrimal(self.phys_pb, 0.0)

    def update(self, convManager):
        """Update the physical state.

        Arguments:
            convManager (ConvergenceManager): Object that manages the
                convergency criteria.
        """

        if self._get("CONTACT", "ALGO_RESO_GEOM") == "POINT_FIXE":
            geom_diff = self.phys_state.primal_step - self.geom_step
            self.geom_step = self.phys_state.primal_step.duplicate()
        else:
            geom_diff = self.phys_state.createPrimal(self.phys_pb, 0.0)

        convManager.evalGeometricResidual(geom_diff)

        self.geom = (
            self.phys_pb.getMesh().getCoordinates()
            + self.phys_state.primal
            + self.phys_state.primal_step
        )

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
    def solve(self):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        convManager = self.get_feature(SOP.ConvergenceManager)
        iter_geom = convManager.setdefault("ITER_GEOM")
        iter_geom.value = -1
        logManager = self.createLoggingManager()
        # logManager.printIntro(self.phys_state.time + self.phys_state.time_step, 1)
        logManager.printConvTableEntries()

        self.geom = self.phys_pb.getMesh().getCoordinates() + self.phys_state.primal

        criteria = self.get_feature(SOP.ConvergenceCriteria)

        while not convManager.isFinished():
            self.current_incr += 1
            iter_geom.value = self.current_incr

            if criteria.has_feature(SOP.Contact):
                criteria.get_feature(SOP.Contact).setPairingCoordinates(self.geom)
            criteria.setLoggingManager(logManager)
            criteria.initialize()

            # Solve current iteration
            self.current_matrix = criteria.solve(self.current_matrix)

            # Update physical state
            self.update(convManager)

            if self._get("CONTACT", "ALGO_RESO_GEOM") == "POINT_FIXE":
                logManager.printConvTableRow(
                    [self.current_incr, " ", " ", convManager.get("RESI_GEOM"), "POINT_FIXE"]
                )

        if not convManager.isConverged():
            raise ConvergenceError("MECANONLINE9_9")

        deleteTemporaryObjects()

        logManager.printConvTableEnd()

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        args = self.param
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword])[0].get(parameter, default)

        return args.get(keyword, default)
