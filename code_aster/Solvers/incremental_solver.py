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

from ..Objects import AssemblyMatrixDisplacementReal, DiscreteComputation
from ..Supervis import IntegrationError
from ..Utilities import no_new_attributes, profile
from .base_features import EventSource
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP


class Residuals:
    """Container to store intermediate field.

    Attributes:
        resi (FieldOnNodesReal): Global residual.
        resi_int (FieldOnNodesReal):  Internal residual.
        resi_ext (FieldOnNodesReal): External residual.
        resi_dual (FieldOnNodesReal): Dirichlet reactions.
        resi_stress (FieldOnNodesReal): Internal forces.
        resi_cont (FieldOnNodesReal): Contact residual.
    """

    resi = resi_int = resi_ext = resi_dual = resi_stress = resi_cont = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def update(self):
        if self.resi:
            self.resi.updateValuePointers()
        if self.resi_int:
            self.resi_int.updateValuePointers()
        if self.resi_ext:
            self.resi_ext.updateValuePointers()
        if self.resi_dual:
            self.resi_dual.updateValuePointers()
        if self.resi_stress:
            self.resi_stress.updateValuePointers()
        if self.resi_cont:
            self.resi_cont.updateValuePointers()


class IncrementalSolver(SolverFeature, EventSource):
    """Solve an iteration."""

    provide = SOP.IncrementalSolver | SOP.EventSource
    required_features = [SOP.PhysicalProblem, SOP.PhysicalState, SOP.LinearSolver, SOP.LineSearch]
    optional_features = [SOP.Contact, SOP.ConvergenceManager]

    _data = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()
        self._data = {}

    def notifyObservers(self, convManager, matrix_type):
        """Notify all observers about the convergence.

        Arguments:
            convManager (ConvergenceManager): Object that holds the criteria values.
            matrix_type (str): Type of matrix used.
        """
        self._data = convManager.getParameters()
        self._data["matrix"] = matrix_type
        self._data["isConverged"] = convManager.isConverged()
        super().notifyObservers()

    def get_state(self):
        """Returns the current residuals to be shared with observers."""
        return SOP.IncrementalSolver, self._data

    @profile
    @SolverFeature.check_once
    def computeInternalResidual(self, residuals, scaling=1.0):
        """Compute internal residual R_int(u, Lagr).

            R_int(u, Lagr) = [B^t.Sig(u) + B^t.Lagr, B^t.primal-primal_impo]

        Arguments:
            residuals (Residuals):
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0)

        Returns:
            FieldOnCellsReal: internal state variables (VARI_ELGA)
            FieldOnCellsReal: Cauchy stress tensor (SIEF_ELGA)
        """
        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute internal forces (B^t.stress)
        _, codret, internVar, stress, r_stress = disc_comp.getInternalForces(
            self.phys_state.primal,
            self.phys_state.primal_step,
            self.phys_state.stress,
            self.phys_state.internVar,
            self.phys_state.time,
            self.phys_state.time_step,
            self.phys_state.getState(-1).externVar,
            self.phys_state.externVar,
        )

        residuals.resi_stress = r_stress
        r_int = r_stress

        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        if self.phys_pb.getDOFNumbering().useLagrangeDOF():
            primal_curr = self.phys_state.primal + self.phys_state.primal_step

            # Compute kinematic forces (B^t.Lagr_curr)
            dualizedBC_forces = disc_comp.getDualForces(primal_curr)
            residuals.resi_dual = dualizedBC_forces

            # Compute dualized BC (B^t.primal_curr - primal_impo)
            # Compute dualized BC (B^t.primal_curr)
            dualizedBC_disp = disc_comp.getDualDisplacement(primal_curr, scaling)

            # Imposed dualized BC (primal_impo)
            time_curr = self.phys_state.time + self.phys_state.time_step
            dualizedBC_impo = disc_comp.getImposedDualBC(time_curr)

            r_int += dualizedBC_forces + dualizedBC_disp - dualizedBC_impo
        else:
            residuals.resi_dual = self.phys_state.createPrimal(self.phys_pb, 0.0)

        residuals.resi_int = r_int

        return internVar, stress

    @profile
    @SolverFeature.check_once
    def computeExternalResidual(self, residuals):
        """Compute external residual R_ext(u, Lagr)

            R_ext(u, Lagr) = [(f(u),v), 0]

        Arguments:
            residuals (Residuals):

        Returns:
            FieldOnNodesReal: External residual.
        """
        # Main object for discrete computation
        disc_comp = DiscreteComputation(self.phys_pb)

        # Compute Neumann forces
        neumann_forces = disc_comp.getNeumannForces(
            self.phys_state.time + self.phys_state.time_step, varc_curr=self.phys_state.externVar
        )
        volum_forces = disc_comp.getVolumetricForces(
            self.phys_state.time + self.phys_state.time_step, varc_curr=self.phys_state.externVar
        )
        residuals.resi_ext = neumann_forces + volum_forces

    @profile
    @SolverFeature.check_once
    def computeContactResidual(self, residuals):
        """Compute contact residual R_cont(u, Lagr)

        Arguments:
            residuals (Residuals):
        """
        contact_manager = self.get_feature(SOP.Contact, optional=True)
        if contact_manager:
            disc_comp = DiscreteComputation(self.phys_pb)

            # Compute contact forces
            contact_forces = disc_comp.getContactForces(
                contact_manager.getPairingCoordinates(),
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.phys_state.time,
                self.phys_state.time_step,
                contact_manager.data(),
                contact_manager.coef_cont,
                contact_manager.coef_frot,
            )
        else:
            contact_forces = self.phys_state.createPrimal(self.phys_pb, 0.0)

        residuals.resi_cont = contact_forces

    @profile
    @SolverFeature.check_once
    def computeResidual(self, scaling=1.0):
        """Compute R(u, Lagr) = - (Rint(u, Lagr) + Rcont(u, Lagr) - Rext(u, Lagr)).

        This is not the true residual but the opposite.

        Arguments:
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0).

        Returns:
            tuple(Residuals, FieldOnCellsReal, FieldOnCellsReal):
            Tuple with residuals, internal state variables (VARI_ELGA),
            Cauchy stress tensor (SIEF_ELGA).
        """
        residuals = Residuals()
        internVar, stress = self.computeInternalResidual(residuals, scaling)
        self.computeExternalResidual(residuals)
        self.computeContactResidual(residuals)
        # Compute residual
        residuals.resi = -(residuals.resi_int + residuals.resi_cont - residuals.resi_ext)

        return residuals, internVar, stress

    @profile
    @SolverFeature.check_once
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
            time_curr = self.phys_state.time + self.phys_state.time_step
            matr_elem_rigi = disc_comp.getLinearStiffnessMatrix(time=time_curr, with_dual=False)
            codret = 0
        elif matrix_type == "PRED_TANGENTE":
            _, codret, matr_elem_rigi = disc_comp.getPredictionTangentStiffnessMatrix(
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.phys_state.stress,
                self.phys_state.internVar,
                self.phys_state.time,
                self.phys_state.time_step,
                self.phys_state.getState(-1).externVar,
                self.phys_state.externVar,
            )
        elif matrix_type == "TANGENTE":
            _, codret, matr_elem_rigi = disc_comp.getTangentStiffnessMatrix(
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.phys_state.stress,
                self.phys_state.internVar,
                self.phys_state.time,
                self.phys_state.time_step,
                self.phys_state.getState(-1).externVar,
                self.phys_state.externVar,
            )
        else:
            raise RuntimeError("Matrix not supported: %s" % (matrix_type))

        # Compute dual matrix
        matr_elem_dual = disc_comp.getDualStiffnessMatrix()

        return codret, matr_elem_rigi, matr_elem_dual

    @profile
    @SolverFeature.check_once
    def computeContactJacobian(self):
        """Compute K(u) = d(Rcont(u) ) / du

        Returns:
           ElementaryMatrixDisplacementReal: Contact matrix.
        """
        contact_manager = self.get_feature(SOP.Contact, optional=True)
        if contact_manager:
            # Main object for discrete computation
            disc_comp = DiscreteComputation(self.phys_pb)

            matr_elem_cont = disc_comp.getContactMatrix(
                contact_manager.getPairingCoordinates(),
                self.phys_state.primal,
                self.phys_state.primal_step,
                self.phys_state.time,
                self.phys_state.time_step,
                contact_manager.data(),
                contact_manager.coef_cont,
                contact_manager.coef_frot,
            )

            return matr_elem_cont

        return None

    @profile
    @SolverFeature.check_once
    def computeJacobian(self, matrix_type):
        """Compute K(u) = d(Rint(u) - Rext(u)) / du

        Arguments:
            matrix_type (str): type of matrix used.

        Returns:
            AssemblyMatrixDisplacementReal: Jacobian matrix.
        """
        # Compute elementary matrix
        codret, matr_elem_rigi, matr_elem_dual = self.computeInternalJacobian(matrix_type)
        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        matr_elem_cont = self.computeContactJacobian()

        # Assemble matrix
        jacobian = AssemblyMatrixDisplacementReal(self.phys_pb)

        jacobian.addElementaryMatrix(matr_elem_rigi)
        jacobian.addElementaryMatrix(matr_elem_dual)
        jacobian.addElementaryMatrix(matr_elem_cont)

        jacobian.assemble()

        return jacobian

    @profile
    @SolverFeature.check_once
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
        residuals, internVar, stress = self.computeResidual(scaling)

        # evaluate convergence
        convManager = self.get_feature(SOP.ConvergenceManager)
        convManager.evalNormResidual(residuals)
        self.notifyObservers(convManager, matrix_type)

        if not convManager.isConverged():
            # Time at end of current step
            time_curr = self.phys_state.time + self.phys_state.time_step

            # Compute Dirichlet BC:
            disc_comp = DiscreteComputation(self.phys_pb)
            primal_curr = self.phys_state.primal + self.phys_state.primal_step
            diriBCs = disc_comp.getIncrementalDirichletBC(time_curr, primal_curr)

            # Solve linear system
            linear_solver = self.get_feature(SOP.LinearSolver)
            if not stiffness.isFactorized():
                linear_solver.factorize(stiffness, raiseException=True)
            solution = linear_solver.solve(residuals.resi, diriBCs)

            # Use line search
            lineSearch = self.get_feature(SOP.LineSearch)

            primal_incr = lineSearch.solve(solution)
        else:
            primal_incr = self.phys_state.createPrimal(self.phys_pb, 0.0)

        return primal_incr, internVar, stress, stiffness
