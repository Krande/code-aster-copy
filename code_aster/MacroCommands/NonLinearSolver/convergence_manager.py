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

from ...NonLinear import SolverFeature
from ...NonLinear import SolverOptions as SOP
from ...Objects import DiscreteComputation
from ...Utilities import MPI, no_new_attributes, profile


class ConvergenceManager(SolverFeature):
    """Object that decides about the convergence status."""

    provide = SOP.ConvergenceManager
    required_features = [SOP.PhysicalProblem, SOP.PhysicalState]

    criteria = values = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()
        self.criteria = {}
        self.values = {}

    def initialize(self):
        """Initialize the object for a new iteration."""
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
    @SolverFeature.check_once
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
    @SolverFeature.check_once
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
    @SolverFeature.check_once
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
    @SolverFeature.check_once
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

        if not self.values:
            if not self.criteria:
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
