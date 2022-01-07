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

from ...Utilities import no_new_attributes, profile, MPI
from ...Objects import DiscreteComputation


class ConvergenceManager:
    """Object that decides about the convergence status.

    Arguments:
        epsilon (float): Expected precision.
        test_type (str): Type of criteria (absolute or relative).
        phys_pb (PhysicalProblem): Physical problem.
        phys_state (PhysicalState): Physical state.
    """
    converged = epsilon = residual = test_type = None
    phys_state = phys_pb = elem_comp = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, epsilon, test_type, phys_pb, phys_state):
        # store convergence parameters : RESI_GLOB_RELA, RESI_GLOB_MAXI...
        self.converged = False
        self.epsilon = epsilon
        self.phys_pb = phys_pb
        self.phys_state = phys_state
        assert test_type in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI"), test_type
        self.test_type = test_type

    @profile
    def getDirichletResidual(self, residual):
        """Return the residual with Dirichlet imposed values.

        Arguments:
            residual (FieldOnNodesReal): Residual.

        Returns:
            FieldOnNodesReal: Residual changed in place.
        """
        dofNume = self.phys_pb.getDOFNumbering()

        # maybe not really efficient
        if dofNume.hasDirichletBC():
            time_curr = self.phys_state.time + self.phys_state.time_step
            displ_curr = self.phys_state.displ + self.phys_state.displ_incr
            disc_comp = DiscreteComputation(self.phys_pb)
            diriBCs = disc_comp.incrementalDirichletBC(time_curr, displ_curr)
            eliminatedDofs = dofNume.getDirichletBCDOFs()
            nbElimination = len(eliminatedDofs)
            assert residual.size() == nbElimination

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

        dofNume = self.phys_pb.getDOFNumbering()

        eliminatedDofs = dofNume.getDirichletBCDOFs()
        nb_dofs = len(eliminatedDofs)

        for ieq in range(nb_dofs):
            f_int = 0.0
            f_ext = 0.0
            f_varc = 0.0
            if eliminatedDofs[ieq] == 1:
                f_int = - residuals.resi_int[ieq]
            else:
                f_int = residuals.resi_dual[ieq]
                f_ext = residuals.resi_ext[ieq]

            value = abs(f_int - f_ext) + abs(f_varc)

            if scaling < value:
                scaling = value

        return MPI.COMM_WORLD.allreduce(scaling, MPI.MAX)

    @profile
    def getNormResidual(self, residuals):
        """Returns a dictionnary of criteria

        Arguments:
            residuals (ResiState): Collections of residuals.

        Returns:
            dict: convergence criteria (maxi, relative).
        """

        criteria = {}

        residual = self.getDirichletResidual(residuals.resi)

        criteria["RESI_GLOB_MAXI"] = residual.norm("NORM_INFINITY")

        scaling = self.getRelativeScaling(residuals)

        if scaling == 0.0:
            criteria["RESI_GLOB_RELA"] = -1.0
        else:
            criteria["RESI_GLOB_RELA"] = criteria["RESI_GLOB_MAXI"] / scaling

        return criteria

    def isConverged(self, residuals):
        """Tell if the *residual* is acceptable.

        Arguments:
            residual (ResiState): Collection of residuals.

        Returns:
            bool: *True* if converged, *False* otherwise.
        """
        self.residual = self.getNormResidual(residuals)
        test_type = self.test_type
        if test_type == "RESI_GLOB_RELA" and self.residual[test_type] < -0.5:
            test_type = "RESI_GLOB_MAXI"
        self.converged = self.residual[test_type] < self.epsilon
        return self.converged

    def hasConverged(self):
        """Tell if the last iteration is converged.

        Returns:
            bool: *True* if converged, *False* otherwise.
        """
        return self.converged
