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

from ...Objects import DiscreteComputation, AssemblyMatrixDisplacementReal
from ...Supervis import IntegrationError
from ...Utilities import no_new_attributes, profile
from .elementary_computation import ElementaryComputation


class ResiState:
    """Container to store intermediate field."""

    resi = resi_int = resi_ext = resi_dual = resi_stress = resi_cont = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def update(self):
        if self.resi is not None:
            self.resi.updateValuePointers()
        if self.resi_int is not None:
            self.resi_int.updateValuePointers()
        if self.resi_ext is not None:
            self.resi_ext.updateValuePointers()
        if self.resi_dual is not None:
            self.resi_dual.updateValuePointers()
        if self.resi_stress is not None:
            self.resi_stress.updateValuePointers()
        if self.resi_cont is not None:
            self.resi_cont.updateValuePointers()


class IncrementalSolver:
    """Solve an iteration."""

    phys_state = phys_pb = convManager = None
    linear_solver = None
    elem_comp = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self.elem_comp = ElementaryComputation()

    def setPhysicalProblem(self, phys_pb):
        """Assign the physical problem.

        Arguments:
            phys_pb (PhysicalProblem): Physical problem
        """
        self.phys_pb = phys_pb
        self.elem_comp.setPhysicalProblem(self.phys_pb)

    def setPhysicalState(self, phys_state):
        """Assign the physical state.

        Arguments:
            phys_state (PhysicalState): Physical state
        """
        self.phys_state = phys_state
        self.elem_comp.setPhysicalState(self.phys_state)

    def setConvergenceCriteria(self, convManager):
        """Set the convergence criteria to be used.

        Arguments:
            convManager (ConvergenceManager): object to be used.
        """
        self.convManager = convManager

    def setLinearSolver(self, solver):
        """Set the linear solver to be used.

        Arguments:
            solver (LinearSolver): a linear solver object
        """
        self.linear_solver = solver

    @profile
    def computeInternalResidual(self, timeFieldEndStep, scaling=1.0):
        """Compute internal residual R_int(u, Lagr).

            R_int(u, Lagr) = [B^t.Sig(u) + B^t.Lagr, B^t.displ-displ_impo]

        Arguments:
            timeFieldEndStep (ConstantFieldOnCellsReal): field for time at end of time step
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0)

        Returns:
            ResiState: residuals (stress, dual, internal)
            FieldOnCellsReal: internal state variables (VARI_ELGA)
            FieldOnCellsReal: Cauchy stress tensor (SIEF_ELGA)
        """
        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute internal forces (B^t.stress)
        codret, internVar, stress, r_stress = self.elem_comp.computeInternalForces(
            self.phys_state.time_field, timeFieldEndStep
        )

        resi_state = ResiState()
        resi_state.resi_stress = r_stress
        r_int = r_stress

        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        # Update displacement
        timeEndStep = self.phys_state.time + self.phys_state.time_step
        displ_curr = self.phys_state.displ + self.phys_state.displ_incr

        if self.phys_pb.getDOFNumbering().useLagrangeMultipliers():
            # Compute kinematic forces (B^t.Lagr_curr)
            dualizedBC_forces = disc_comp.dualReaction(displ_curr)
            resi_state.resi_dual = dualizedBC_forces

            # Compute dualiazed BC (B^t.displ_curr - displ_impo)
            # Compute dualiazed BC (B^t.displ_curr)
            dualizedBC_disp = disc_comp.dualDisplacement(displ_curr, scaling)

            # Imposed dualisez BC (displ_impo)
            dualizedBC_impo = disc_comp.imposedDisplacement(timeEndStep)

            r_int += dualizedBC_forces + dualizedBC_disp - dualizedBC_impo
        else:
            resi_state.resi_dual = self.phys_state.createDisplacement(self.phys_pb, 0.0)

        resi_state.resi_int = r_int

        return resi_state, internVar, stress

    @profile
    def computeExternalResidual(self):
        """Compute external residual R_ext(u, Lagr)

            R_ext(u, Lagr) = [(f(u),v), 0]

        Returns:
            FieldOnNodesReal: External residual.
        """
        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute neuamnn forces
        times = [self.phys_state.time, self.phys_state.time_step, 0.0]
        neumann_forces = disc_comp.neumann(times)
        r_ext = neumann_forces

        return r_ext

    @profile
    def computeResidual(self, timeFieldEndStep, scaling=1.0):
        """Compute R(u, Lagr) = - (Rint(u, Lagr) - Rext(u, Lagr)).

        This is not the true residual but the opposite.

        Arguments:
            timeFieldEndStep (ConstantFieldOnCellsReal): field for time at end of time step
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0).

        Returns:
            tuple (ResiState(), FieldOnCellsReal, FieldOnCellsReal):
            Tuple with residuals, internal state variables (VARI_ELGA),
            Cauchy stress tensor (SIEF_ELGA).
        """

        # Compute internal residual
        resi_state, internVar, stress = self.computeInternalResidual(timeFieldEndStep, scaling)

        # Compute external residual
        resi_state.resi_ext = self.computeExternalResidual()

        # Compute residual
        resi_state.resi = -(resi_state.resi_int - resi_state.resi_ext)

        # clean temporary memory - too many objects are not destroyed in fortran
        deleteTemporaryObjects()

        return resi_state, internVar, stress

    @profile
    def computeJacobian(self, matrix_type, timeFieldEndStep):
        """Compute K(u) = d(Rint(u) - Rext(u)) / du

        Arguments:
            matrix_type (str): type of matrix used.
            timeFieldEndStep (ConstantFieldOnCellsReal): field for time at end of time step

        Returns:
            AssemblyMatrixDisplacementReal: Jacobian matrix.
        """
        # Compute elementary matrix
        matr_elem_diri = None
        if matrix_type in ("PRED_ELASTIQUE", "ELASTIQUE"):
            codret, matr_elem = self.elem_comp.computeElasticStiffnessMatrix()
        elif matrix_type == "PRED_TANGENTE":
            codret, matr_elem = self.elem_comp.computeTangentPredictionMatrix(
                self.phys_state.time_field, timeFieldEndStep
            )
            matr_elem_diri = self.elem_comp.computeDualStiffnessMatrix()
        elif matrix_type == "TANGENTE":
            codret, matr_elem = self.elem_comp.computeTangentStiffnessMatrix(
                self.phys_state.time_field, timeFieldEndStep
            )
            matr_elem_diri = self.elem_comp.computeDualStiffnessMatrix()
        else:
            raise RuntimeError("Matrix not supported: %s" % (matrix_type))

        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        # Assembly matrix
        jacobian = AssemblyMatrixDisplacementReal(self.phys_pb)
        jacobian.addElementaryMatrix(matr_elem)
        if matr_elem_diri:
            if matr_elem_diri.hasElementaryTerms():
                jacobian.addElementaryMatrix(matr_elem_diri)
        jacobian.assemble()

        # clean temporary memory - too many objects are not destroyed in fortran
        deleteTemporaryObjects()

        return jacobian

    @profile
    def lineSearch(self, field):
        """Apply linear search.

        Arguments:
            field (FieldOnNodes): DIsplacement field.

        Returns:
            FieldOnNodes: Accelerated field by linear search.
        """
        return 1.0 * field

    @profile
    def solve(self, matrix_type, timeFieldEndStep, matrix=None):
        """Solve the iteration.

        Arguments:
            matrix_type (str): type of matrix used.
            matrix (AssemblyMatrixDisplacementReal, optional): Stiffness matrix
                to be reused.
            timeFieldEndStep (ConstantFieldOnCellsReal): field for time at end of time step

        Returns:
            tuple (FieldOnNodesReal, FieldOnCellsReal, FieldOnCellsReal,
            AssemblyMatrixDisplacementReal): Tuple with incremental displacement
            (DEPL), internal state variables (VARI_ELGA), Cauchy stress (SIEF_ELGA),
            Jacobian matrix used (if computed).
        """

        # we need the matrix to have scaling factor for Lagrange
        # On peut sûrement faire une optimisation ici quand la matrice tangente
        # est calculée (donc le comportement). Il faudrait dans ce cas et
        # uniquement dans ce cas calculer les forces internes avec RAPH_MECA
        # et pas tout recalculer (a voir plus tard)
        if not matrix:
            stiffness = self.computeJacobian(matrix_type, timeFieldEndStep)
        else:
            stiffness = matrix

        # Get scaling for Lagrange (j'aimerais bien m'en débarasser de là)
        scaling = stiffness.getLagrangeScaling()

        # compute residual
        resi_state, internVar, stress = self.computeResidual(timeFieldEndStep, scaling)

        # clean temporary memory - too many objects are not destroyed in fortran
        deleteTemporaryObjects()

        # Time at end of current step
        timeEndStep = self.phys_state.time + self.phys_state.time_step

        if not self.convManager.isConverged(resi_state):
            # compute Dirichlet BC:
            disc_comp = DiscreteComputation(self.phys_pb)
            displ_curr = self.phys_state.displ + self.phys_state.displ_incr
            diriBCs = disc_comp.incrementalDirichletBC(timeEndStep, displ_curr)

            # solve linear system
            if not stiffness.isFactorized():
                self.linear_solver.factorize(stiffness)
            solution = self.linear_solver.solve(resi_state.resi, diriBCs)

            # use line search
            displ_incr = self.lineSearch(solution)
        else:
            displ_incr = self.phys_state.createDisplacement(self.phys_pb, 0.0)

        # clean temporary memory - too many objects are not destroyed in fortran
        deleteTemporaryObjects()

        return displ_incr, internVar, stress, stiffness
