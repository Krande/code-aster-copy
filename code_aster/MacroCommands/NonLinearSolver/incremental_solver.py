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
from .contact_solver import ContactSolver

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
    linear_solver = contact_solver = None
    __setattr__ = no_new_attributes(object.__setattr__)

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

    def setContactSolver(self, solver):
        """Set the contact solver to be used.

        Arguments:
            solver (ContactSolver): a contact solver object
        """
        self.contact_solver = solver

    @profile
    def computeInternalResidual(self, scaling=1.0):
        """Compute internal residual R_int(u, Lagr).

            R_int(u, Lagr) = [B^t.Sig(u) + B^t.Lagr, B^t.primal-primal_impo]

        Arguments:
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0)

        Returns:
            ResiState: residuals (stress, dual, internal)
            FieldOnCellsReal: internal state variables (VARI_ELGA)
            FieldOnCellsReal: Cauchy stress tensor (SIEF_ELGA)
        """
        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute internal forces (B^t.stress)
        _, codret, internVar, stress, r_stress = disc_comp.computeInternalForces(
            self.phys_state.primal,
            self.phys_state.primal_step,
            self.phys_state.stress,
            self.phys_state.internVar,
            self.phys_state.time,
            self.phys_state.time_step
        )

        resi_state = ResiState()
        resi_state.resi_stress = r_stress
        r_int = r_stress

        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        if self.phys_pb.getDOFNumbering().useLagrangeMultipliers():
            primal_curr = self.phys_state.primal + self.phys_state.primal_step

            # Compute kinematic forces (B^t.Lagr_curr)
            dualizedBC_forces = disc_comp.dualReaction(primal_curr)
            resi_state.resi_dual = dualizedBC_forces

            # Compute dualiazed BC (B^t.primal_curr - primal_impo)
            # Compute dualiazed BC (B^t.primal_curr)
            dualizedBC_disp = disc_comp.dualDisplacement(primal_curr, scaling)

            # Imposed dualisez BC (primal_impo)
            time_curr = self.phys_state.time + self.phys_state.time_step
            dualizedBC_impo = disc_comp.imposedDualBC(time_curr)

            r_int += dualizedBC_forces + dualizedBC_disp - dualizedBC_impo
        else:
            resi_state.resi_dual = self.phys_state.createPrimal(
                self.phys_pb, 0.0)

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
        neumann_forces = disc_comp.neumann(self.phys_state.time + self.phys_state.time_step,
                                            0.0, 0.0)

        return neumann_forces

    @profile
    def computeContactResidual(self):
        """Compute contact residual R_cont(u, Lagr)

        Returns:
            FieldOnNodesReal: contact residual.
        """

        if self.contact_solver.enable:
            disc_comp = DiscreteComputation(self.phys_pb)

            # Compute contact forces
            contact_forces = disc_comp.contactForces(
                self.contact_solver.getPairingCoordinates(),
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.contact_solver.data())
        else:
            contact_forces = self.phys_state.createPrimal(
                self.phys_pb, 0.0)

        return contact_forces

    @profile
    def computeResidual(self, scaling=1.0):
        """Compute R(u, Lagr) = - (Rint(u, Lagr) + Rcont(u, Lagr) - Rext(u, Lagr)).

        This is not the true residual but the opposite.

        Arguments:
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0).

        Returns:
            tuple (ResiState(), FieldOnCellsReal, FieldOnCellsReal):
            Tuple with residuals, internal state variables (VARI_ELGA),
            Cauchy stress tensor (SIEF_ELGA).
        """

        # Compute internal residual
        resi_state, internVar, stress = self.computeInternalResidual(scaling)

        # Compute external residual
        resi_state.resi_ext = self.computeExternalResidual()

        # Compute contact residual
        resi_state.resi_cont = self.computeContactResidual()

        # Compute residual
        resi_state.resi = -(resi_state.resi_int +
                            resi_state.resi_cont - resi_state.resi_ext)

        # clean temporary memory - too many objects are not destroyed in fortran
        deleteTemporaryObjects()

        return resi_state, internVar, stress

    @profile
    def computeInternalJacobian(self, matrix_type):
        """Compute K(u) = d(Rint(u)) / du

        Arguments:
            matrix_type (str): type of matrix used.

        Returns:
            int : error code flag
            ElementaryMatrixDisplacementReal: rigidity matrix.
            ElementaryMatrixDisplacementReal: dual matrix.
        """
        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute rigidity matrix
        if matrix_type in ("PRED_ELASTIQUE", "ELASTIQUE"):
            matr_elem_rigi = disc_comp.elasticStiffnessMatrix()
            codret = 0
        elif matrix_type == "PRED_TANGENTE":
            _, codret, matr_elem_rigi = disc_comp.computeTangentPredictionMatrix(
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.phys_state.stress,
                self.phys_state.internVar,
                self.phys_state.time,
                self.phys_state.time_step
            )
        elif matrix_type == "TANGENTE":
            _, codret, matr_elem_rigi = disc_comp.computeTangentStiffnessMatrix(
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.phys_state.stress,
                self.phys_state.internVar,
                self.phys_state.time,
                self.phys_state.time_step
            )
        else:
            raise RuntimeError("Matrix not supported: %s" % (matrix_type))

        # Compute dual matrix
        matr_elem_dual = disc_comp.dualStiffnessMatrix()

        return codret, matr_elem_rigi, matr_elem_dual

    @profile
    def computeContactJacobian(self):
        """Compute K(u) = d(Rcont(u) ) / du

        Returns:
           ElementaryMatrixDisplacementReal: Contact matrix.
        """

        if self.contact_solver.enable:
            # Main object for discrete computation
            disc_comp = DiscreteComputation(self.phys_pb)

            matr_elem_cont = disc_comp.contactMatrix(
                self.contact_solver.getPairingCoordinates(),
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.contact_solver.data()
            )

            return matr_elem_cont

        return None

    @profile
    def computeJacobian(self, matrix_type):
        """Compute K(u) = d(Rint(u) - Rext(u)) / du

        Arguments:
            matrix_type (str): type of matrix used.

        Returns:
            AssemblyMatrixDisplacementReal: Jacobian matrix.
        """
        # Compute elementary matrix

        codret, matr_elem_rigi, matr_elem_dual = self.computeInternalJacobian(
            matrix_type)

        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        matr_elem_cont = self.computeContactJacobian()

        # Assemble matrix
        jacobian = AssemblyMatrixDisplacementReal(self.phys_pb)

        jacobian.addElementaryMatrix(matr_elem_rigi)
        jacobian.addElementaryMatrix(matr_elem_dual)
        jacobian.addElementaryMatrix(matr_elem_cont)

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
    def solve(self, matrix_type, matrix=None):
        """Solve the iteration.

        Arguments:
            matrix_type (str): type of matrix used.
            matrix (AssemblyMatrixDisplacementReal, optional): Stiffness matrix
                to be reused.

        Returns:
            tuple (FieldOnNodesReal, FieldOnCellsReal, FieldOnCellsReal,
            AssemblyMatrixDisplacementReal): Tuple with incremental primal,
            internal state variables (VARI_ELGA), Cauchy stress (SIEF_ELGA),
            Jacobian matrix used (if computed).
        """

        # we need the matrix to have scaling factor for Lagrange
        if not matrix:
            stiffness = self.computeJacobian(matrix_type)
        else:
            stiffness = matrix

        # Get scaling for Lagrange (j'aimerais bien m'en débarasser de là)
        scaling = stiffness.getLagrangeScaling()

        # compute residual
        resi_state, internVar, stress = self.computeResidual(scaling)

        # evaluate convergence
        self.convManager.evalNormResidual(resi_state)

        if not self.convManager.hasConverged():
            # Time at end of current step
            time_curr = self.phys_state.time + self.phys_state.time_step

            # compute Dirichlet BC:
            disc_comp = DiscreteComputation(self.phys_pb)
            primal_curr = self.phys_state.primal + self.phys_state.primal_step
            diriBCs = disc_comp.incrementalDirichletBC(time_curr, primal_curr)

            # solve linear system
            if not stiffness.isFactorized():
                self.linear_solver.factorize(stiffness)
            solution = self.linear_solver.solve(resi_state.resi, diriBCs)

            # use line search
            primal_incr = self.lineSearch(solution)
        else:
            primal_incr = self.phys_state.createPrimal(self.phys_pb, 0.0)

        # clean temporary memory - too many objects are not destroyed in fortran
        deleteTemporaryObjects()

        return primal_incr, internVar, stress, stiffness