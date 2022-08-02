# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

from ...Supervis import ConvergenceError
from ...Utilities import no_new_attributes, profile
from .convergence_manager import ConvergenceManager
from .geometric_solver import GeometricSolver
from .logging_manager import LoggingManager


class StepSolver:
    """Solves a step, loops on iterations.
    """

    current_incr = phys_state = None
    phys_pb = linear_solver = None
    param = None
    geom = geom_step = current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self.current_incr = 0

    def setParameters(self, param):
        """Set parameters from user keywords.

        Arguments:
            param (dict) : user keywords.
        """
        self.param = param

    def setPhysicalProblem(self, phys_pb):
        """Assign the physical problem.

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
        """
        self.phys_pb = phys_pb

    def setPhysicalState(self, phys_state):
        """Assign the physical state.

        Arguments:
            phys_state (PhysicalState): Physical state
        """
        self.phys_state = phys_state
        self.phys_state.primal_step = self.phys_state.createPrimal(self.phys_pb, 0.0)
        self.geom_step = self.phys_state.createPrimal(self.phys_pb, 0.0)

    def getPhysicalState(self):
        """Get the physical state.

        Returns:
            PhysicalState: Currently used Physical state
        """
        return self.phys_state

    def setLinearSolver(self, linear_solver):
        """Set the linear solver to be used.

        Arguments:
            linear_solver (LinearSolver): Linear solver object.
        """
        self.linear_solver = linear_solver

    def update(self, convManager):
        """Update the physical state.

        Arguments:
            primal_step (FieldOnNodes): Displacement step.
            internVar (FieldOnCells): Internal state variables.
            sigma (FieldOnCells): Stress field.
            convManager (ConvergenceManager): Object that manages the
                convergency criteria.
        """

        if self._get("CONTACT", "ALGO_RESO_GEOM") == "POINT_FIXE":
            geom_diff = self.phys_state.primal_step - self.geom_step
            self.geom_step = self.phys_state.primal_step.duplicate()
        else:
            geom_diff = self.phys_state.createPrimal(self.phys_pb, 0.0)

        convManager.evalGeometricResidual(geom_diff)

        self.geom = self.phys_pb.getMesh().getCoordinates() \
            + self.phys_state.primal + self.phys_state.primal_step

    def createConvergenceManager(self):
        """Return an object that holds convergence criteria.

        Returns:
            ConvergenceManager: object that manages the convergency criteria.
        """

        convMana = ConvergenceManager(self.phys_pb, self.phys_state)

        convMana.addCriteria("RESI_GEOM", self._get("CONTACT", "RESI_GEOM", 10e150))

        return convMana

    def createGeometricSolver(self):
        """Return a solver for the next geometric iteration.

        Returns:
            GeometricSolver: object to solve the next iteration.
        """

        geom_solver = GeometricSolver()
        geom_solver.setParameters(self.param)
        geom_solver.setPhysicalProblem(self.phys_pb)
        geom_solver.setPhysicalState(self.phys_state)
        geom_solver.setLinearSolver(self.linear_solver)
        geom_solver.setCoordinates(self.geom)

        return geom_solver

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

    def hasFinished(self):
        """Tell if there are iterations to be computed.

        Returns:
            bool: *True* if there is no iteration to be computed, *False* otherwise.
        """

        if self.current_incr == 0:
            return False

        reac_geom = self._get("CONTACT", "REAC_GEOM")

        if reac_geom == "AUTOMATIQUE":
            nb_iter = self._get("CONTACT", "ITER_GEOM_MAXI")
        elif reac_geom == "CONTROLE":
            nb_iter = self._get("CONTACT", "NB_ITER_GEOM")
        else:
            nb_iter = 1

        return self.current_incr >= nb_iter

    @profile
    def solve(self):
        """Solve a step.

        Raises:
            *ConvergenceError* exception in case of error.
        """
        convManager = self.createConvergenceManager()
        logManager = self.createLoggingManager()
        logManager.printIntro(self.phys_state.time + self.phys_state.time_step, 1)
        logManager.printConvTableEntries()

        self.geom = self.phys_pb.getMesh().getCoordinates() + self.phys_state.primal

        while not self.hasFinished() and not convManager.hasConverged():
            geometric = self.createGeometricSolver()
            geometric.setLoggingManager(logManager)

            # Solve current iteration
            self.current_matrix = geometric.solve(self.current_matrix)

            # Update physical state
            self.update(convManager)

            self.current_incr += 1

            if self._get("CONTACT", "ALGO_RESO_GEOM") == "POINT_FIXE":
                logManager.printConvTableRow(
                    [
                        self.current_incr, " ", " ",
                        convManager.getCriteria("RESI_GEOM"), "POINT_FIXE",
                    ]
                )

        if not convManager.hasConverged():
            raise ConvergenceError("MECANONLINE9_9")

        deleteTemporaryObjects()

        logManager.printConvTableEnd()

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        if parameter is not None:
            if keyword in self.param and self.param.get(keyword) is not None:
                return self.param.get(keyword).get(parameter, default)
            else:
                return default

        return self.param.get(keyword, default)
