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

"""
:py:class:`DiscreteComputation` --- DiscreteComputation object
******************************************************
"""

from libaster import (
    DiscreteComputation,
    AssemblyMatrixDisplacementReal,
    AssemblyMatrixTemperatureReal,
    AssemblyMatrixPressureComplex,
)

from ..Utilities import injector, profile


@injector(DiscreteComputation)
class ExtendedDiscreteComputation:
    def __getinitargs__(self):
        """Returns the argument required to reinitialize a MaterialField
        object during unpickling.
        """
        return (self.getPhysicalProblem(),)

    @profile
    def getDirichletBC(self, time=0.0):
        """Return the imposed displacement vector used to remove imposed DDL

        Arguments:
              time (float): Current time (default 0.0)

        Returns:
              FieldOnNodes: imposed BC vector
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            return self.getMechanicalDirichletBC(time)
        elif model.isThermal():
            return self.getThermalDirichletBC(time)
        elif model.isAcoustic():
            return self.getAcousticDirichletBC(time)
        else:
            raise RuntimeError("Unknown physics")

    @profile
    def getLinearStiffnessMatrix(
        self, time=0.0, fourierMode=-1, groupOfCells=[], with_dual=True, assembly=False
    ):
        """Return the elementary matrices for stiffness matrix depending of the physic.
        Option RIGI_MECA or RIGI_THER or RIGI_ACOU.

        Arguments:
              time (float): Current time for external state variable evaluation.
                Only needed if material depends on time (default: 0.0)
              fourierMode (int): Fourier mode (default: -1)
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it is empty, the full model is used
              with_dual (bool): compute dual terms or not (default: True)
              assembly (bool): assemble elementary matrix (default: False)
        Returns:
              ElementaryMatrix / AssemblyMatrix: (elementary) stiffness matrix depends
                on 'assembly' keyword
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            matr_elem = self.getElasticStiffnessMatrix(time, fourierMode, groupOfCells, with_dual)

            if assembly:
                matr_asse = AssemblyMatrixDisplacementReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isThermal():
            matr_elem = self.getLinearConductivityMatrix(time, fourierMode, groupOfCells, with_dual)

            if assembly:
                matr_asse = AssemblyMatrixTemperatureReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isAcoustic():
            matr_elem = self.getLinearMobilityMatrix(groupOfCells, with_dual)

            if assembly:
                matr_asse = AssemblyMatrixPressureComplex(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse
        else:
            raise RuntimeError("Unknown physic")

        return matr_elem

    @profile
    def getDualStiffnessMatrix(self, assembly=False):
        """Return the elementary matrices for dual stiffness matrix depending of the physic.

        Arguments:
              assembly (bool): assemble elementary matrix (default: False)
        Returns:
              ElementaryMatrix: elementary dual stiffness matrix
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            matr_elem = self.getDualElasticStiffnessMatrix()

            if assembly:
                matr_asse = AssemblyMatrixDisplacementReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isThermal():
            matr_elem = self.getDualLinearConductivityMatrix()

            if assembly:
                matr_asse = AssemblyMatrixTemperatureReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isAcoustic():
            matr_elem = self.getDualLinearMobilityMatrix()

            if assembly:
                matr_asse = AssemblyMatrixPressureComplex(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        else:
            raise RuntimeError("Unknown physic")

        return matr_elem

    @profile
    def getMassMatrix(self, time=0.0, groupOfCells=[], assembly=False):
        """Return the elementary matrices formass matrix depending of the physic.
        Option MASS_MECA or MASS_THER or MASS_ACOU.

        Arguments:
              time (float): Current time for external state variable evaluation.
                Only needed if material depends on time (default: 0.0)
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it is empty, the full model is used
              assembly (bool): assemble elementary matrix (default: False)
        Returns:
              ElementaryMatrix: elementary mass matrix
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            matr_elem = self.getMechanicalMassMatrix(False, time, groupOfCells)

            if assembly:
                matr_asse = AssemblyMatrixDisplacementReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isThermal():
            matr_elem = self.getLinearCapacityMatrix(time, groupOfCells)

            if assembly:
                matr_asse = AssemblyMatrixTemperatureReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isAcoustic():
            matr_elem = self.getCompressibilityMatrix(groupOfCells)

            if assembly:
                matr_asse = AssemblyMatrixPressureComplex(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse
        else:
            raise RuntimeError("Unknown physic")

        return matr_elem
