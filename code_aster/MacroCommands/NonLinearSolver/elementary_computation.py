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

from ...Cata.Language.SyntaxObjects import _F
from ...Commands import ASSE_VECTEUR, CALC_MATR_ELEM, CALCUL, DEFI_LIST_REEL
from ...Utilities import MPI, haveMPI, no_new_attributes, profile
from ...Objects import DiscreteComputation
from libaster import AsterError


class ElementaryComputation:
    """Compute elementary quantities like Dirichlet bc, Neumann loads, ..."""

    phys_state = phys_pb = None

    __setattr__ = no_new_attributes(object.__setattr__)

    def _getLoads(self):
        """Get list of loads (only CHARGE and not Dirichlet) - to delete quickly"""
        names = []
        listLoads = self.phys_pb.getListOfLoads()

        mecar = listLoads.getMechanicalLoadsReal()
        for load in mecar:
            names.append(load)

        mecaf = listLoads.getMechanicalLoadsFunction()
        for load in mecaf:
            names.append(load)

        if haveMPI():
            mecar = listLoads.getParallelMechanicalLoadsReal()
            for load in mecar:
                names.append(load)

            mecaf = listLoads.getParallelMechanicalLoadsFunction()
            for load in mecaf:
                names.append(load)

        return names

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

    @profile
    def computeInternalForces(self, timeFieldBeginStep, timeFieldEndStep):
        """Compute internal forces, stress and internal state variables.

        Arguments:
            timeFieldBeginStep (ConstantFieldOnCellsReal): field time at begin of time step
            timeFieldEndStep (ConstantFieldOnCellsReal): field time at end of time step

        Returns:
            tuple (int, FieldOnCells, FieldOnCells, FieldOnNodes):
            Tuple with 4 objects: exitcode, internal state variables (VARI_ELGA),
            Cauchy stress (SIEF_ELGA), vector of internal forces (`B^T \sigma`).
        """

        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute internal forces
        res = disc_comp.computeInternalForces(
            self.phys_state.displ,
            self.phys_state.displ_incr,
            self.phys_state.stress,
            self.phys_state.internVar,
            timeFieldBeginStep,
            timeFieldEndStep,
        )

        # Check integration of behavior
        codret = MPI.COMM_WORLD.allreduce(res[1], op=MPI.MAX)

        # forces internes assemblees
        assembly_forces = ASSE_VECTEUR(VECT_ELEM=res[4], NUME_DDL=self.phys_pb.getDOFNumbering())

        return codret, res[2], res[3], assembly_forces

    @profile
    def computeElasticStiffnessMatrix(self):
        """Compute elastic matrix (not assembled)

        Returns:
            tuple (int, ElementaryMatrixDisplacementReal): Tuple with exitcode,
            elementary elastic matrix.
        """

        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute elastic stiffness matrix
        elem_matr = disc_comp.elasticStiffnessMatrix()

        # Check integration of behavior
        codret = 0

        return codret, elem_matr

    @profile
    def computeTangentPredictionMatrix(self, timeFieldBeginStep, timeFieldEndStep):
        """Compute tangent prediction matrix (not assembled).

        Arguments:
            timeFieldBeginStep (ConstantFieldOnCellsReal): field time at begin of time step
            timeFieldEndStep (ConstantFieldOnCellsReal): field time at end of time step

        Returns:
            tuple (int, ElementaryMatrixDisplacementReal): Tuple with exitcode,
            elementary tangent prediction matrix.
        """
        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute tangent prediction matrix
        res = disc_comp.computeTangentPredictionMatrix(
            self.phys_state.displ,
            self.phys_state.displ_incr,
            self.phys_state.stress,
            self.phys_state.internVar,
            timeFieldBeginStep,
            timeFieldEndStep,
        )

        # Check integration of behavior
        codret = 0

        return codret, res[2]

    @profile
    def computeTangentStiffnessMatrix(self, timeFieldBeginStep, timeFieldEndStep):
        """Compute tangent matrix (not assembled).

        Arguments:
            timeFieldBeginStep (ConstantFieldOnCellsReal): field time at begin of time step
            timeFieldEndStep (ConstantFieldOnCellsReal): field time at end of time step

        Returns:
            tuple (int, ElementaryMatrixDisplacementReal): Tuple with exitcode,
            elementary tangent matrix.
        """

        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute tangent prediction matrix
        res = disc_comp.computeTangentStiffnessMatrix(
            self.phys_state.displ,
            self.phys_state.displ_incr,
            self.phys_state.stress,
            self.phys_state.internVar,
            timeFieldBeginStep,
            timeFieldEndStep,
        )

        # Check integration of behavior
        codret = MPI.COMM_WORLD.allreduce(res[1], op=MPI.MAX)

        return codret, res[5]

    @profile
    def computeDualStiffnessMatrix(self):
        """Compute dual stiffness matrix for boundary conditions (not assembled).


        Returns:
            ElementaryMatrixDisplacementReal: elementary  matrix.
        """

        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute
        res = disc_comp.dualStiffnessMatrix()

        return res
