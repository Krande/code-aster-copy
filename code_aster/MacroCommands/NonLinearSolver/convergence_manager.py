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

from math import sqrt

from ...Utilities import no_new_attributes, profile, MPI
from ...Objects import DiscreteComputation


class ConvergenceManager:
    """Object that decides about the convergence status.

    Arguments:
        criteria (dict): Expected precision for each criteria.
        values (dict): Current value of criteria.
        phys_pb (PhysicalProblem): Physical problem.
        phys_state (PhysicalState): Physical state.
    """

    criteria = values = None
    phys_state = phys_pb = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, phys_pb, phys_state):
        self.phys_pb = phys_pb
        self.phys_state = phys_state
        self.criteria = {}
        self.values = {}

    def addCriteria(self, criteria, value):
        """Add a convergence criteria to verify

        Arguments:
            criteria (str): name of the criteria.
            value (float): criteria value
        """

        self.criteria[criteria] = value

    def getCriteria(self, criteria):
        """Get current value of the criteria

        Arguments:
            criteria (str): name of the criteria.

        Returns:
            (float): criteria value
        """

        return self.values[criteria]

    @profile
    def getDirichletResidual(self, residual):
        """Return the residual with Dirichlet imposed values.

        Arguments:
            residual (FieldOnNodesReal): Residual.

        Returns:
            FieldOnNodesReal: Residual changed in place.
        """
        loads = self.phys_pb.getListOfLoads()

        # maybe not really efficient
        if loads.hasDirichletBC():
            time_curr = self.phys_state.time + self.phys_state.time_step
            primal_curr = self.phys_state.primal + self.phys_state.primal_step
            disc_comp = DiscreteComputation(self.phys_pb)
            diriBCs = disc_comp.getIncrementalDirichletBC(time_curr, primal_curr)
            eliminatedDofs = self.phys_pb.getDirichletBCDOFs()
            nbElimination = len(eliminatedDofs)
            assert residual.size() == nbElimination

            residual.updateValuePointers()
            diriBCs.updateValuePointers()
            for ieq in range(nbElimination):
                if eliminatedDofs[ieq] == 1:
                    residual[ieq] = diriBCs[ieq]

        return residual

    @profile
    def getRelativeScaling(self, residuals):
        """Returns the scaling fator to compute the relative error

        Arguments:
            residuals (ResiState): Collections of residuals.

        Returns:
            float: scaling factor.
        """

        scaling = 0.0

        eliminatedDofs = self.phys_pb.getDirichletBCDOFs()
        nb_dofs = len(eliminatedDofs)

        residuals.update()

        for ieq in range(nb_dofs):
            f_int = 0.0
            f_ext = 0.0
            f_varc = 0.0
            if eliminatedDofs[ieq] == 1:
                f_int = -residuals.resi_int[ieq]
            else:
                f_int = residuals.resi_dual[ieq]
                f_ext = residuals.resi_ext[ieq]

            value = abs(f_int - f_ext) + abs(f_varc)

            if scaling < value:
                scaling = value

        return MPI.ASTER_COMM_WORLD.allreduce(scaling, MPI.MAX)

    @profile
    def evalNormResidual(self, residuals):
        """Evaluate criteria

        Arguments:
            residuals (ResiState): Collections of residuals.
        """

        residual = self.getDirichletResidual(residuals.resi)

        self.values["RESI_GLOB_MAXI"] = residual.norm("NORM_INFINITY")

        scaling = self.getRelativeScaling(residuals)

        if scaling == 0.0:
            self.values["RESI_GLOB_RELA"] = -1.0
        else:
            self.values["RESI_GLOB_RELA"] = self.values["RESI_GLOB_MAXI"] / scaling

    @profile
    def evalGeometricResidual(self, displ_delta):
        """Evaluate criteria

        Arguments:
            displ_dela (FieldOnNodesReal): variation of displacement.
        """

        # scaling with diagonal of bounding box
        TABG = self.phys_pb.getMesh().getTable("CARA_GEOM")
        x_diag = TABG["X_MAX", 1] - TABG["X_MIN", 1]
        y_diag = TABG["Y_MAX", 1] - TABG["Y_MIN", 1]
        z_diag = TABG["Z_MAX", 1] - TABG["Z_MIN", 1]
        diag = sqrt(pow(x_diag, 2) + pow(y_diag, 2) + pow(z_diag, 2))

        self.values["RESI_GEOM"] = displ_delta.norm("NORM_INFINITY", ["DX", "DY", "DZ"]) / diag

    @profile
    def hasConverged(self):
        """Tell if convergence criteria are verified.

        Returns:
            bool: *True* if converged, *False* otherwise.
        """

        if len(self.values) == 0:
            if len(self.criteria) == 0:
                return True
            else:
                return False

        for crit in self.criteria:
            if crit in self.values:
                if self.values[crit] > self.criteria[crit] or self.values[crit] < 0.0:
                    return False
            else:
                return False

        return True
