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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`DiscreteComputation` --- DiscreteComputation object
******************************************************
"""

from libaster import DiscreteComputation

from ..Utilities import injector, profile


@injector(DiscreteComputation)
class ExtendedDiscreteComputation:

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MaterialField
        object during unpickling.
        """
        return (self.getPhysicalProblem(), )

    @profile
    def getLinearStiffnessMatrix(self, time=0.0, fourierMode=-1, groupOfCells=[], with_dual=True):
        """Return the elementary matrices for stiffness matrix depending of the physic.
            Option RIGI_MECA or RIGI_THER or RIGI_ACOU.

            Arguments:
                  time (float): Current time for external state variable evaluation.
                    Only needed if material depends on time (default: 0.0)
                  fourierMode (int): Fourier mode (default: -1)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it is empty, the full model is used
                  with_dual (bool): compute dual terms or not (default: True)
            Returns:
                  ElementaryMatrix: elementary stiffness matrix
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            return self.getElasticStiffnessMatrix(time, fourierMode, groupOfCells, with_dual)
        elif model.isThermal():
            return self.getLinearConductivityMatrix(time, fourierMode, groupOfCells, with_dual)
        elif model.isAcoustic():
            return self.getLinearMobilityMatrix(groupOfCells, with_dual)
        else:
            raise RuntimeError("Unknown physic")

    @profile
    def getDualStiffnessMatrix(self):
        """Return the elementary matrices for dual stiffness matrix depending of the physic.

            Returns:
                  ElementaryMatrix: elementary dual stiffness matrix
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            return self.getDualElasticStiffnessMatrix()
        elif model.isThermal():
            return self.getDualLinearConductivityMatrix()
        elif model.isAcoustic():
            return self.getDualLinearMobilityMatrix()
        else:
            raise RuntimeError("Unknown physic")

    @profile
    def getMassMatrix(self, time=0.0, groupOfCells=[]):
        """Return the elementary matrices formass matrix depending of the physic.
            Option MASS_MECA or MASS_THER or MASS_ACOU.

            Arguments:
                  time (float): Current time for external state variable evaluation.
                    Only needed if material depends on time (default: 0.0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it is empty, the full model is used
            Returns:
                  ElementaryMatrix: elementary mass matrix
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            return self.getMechanicalMassMatrix(False, time, groupOfCells)
        elif model.isThermal():
            return self.getLinearCapacityMatrix(time, groupOfCells)
        elif model.isAcoustic():
            return self.getCompressibilityMatrix(groupOfCells)
        else:
            raise RuntimeError("Unknown physic")