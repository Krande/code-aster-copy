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

from ..Solvers.residual import Residuals
from ..Supervis import IntegrationError
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
        self,
        time=0.0,
        fourierMode=-1,
        varc_curr=None,
        groupOfCells=[],
        with_dual=True,
        assembly=False,
    ):
        """Return the elementary matrices for stiffness matrix depending of the physic.
        Option RIGI_MECA or RIGI_THER or RIGI_ACOU.

        Arguments:
              time (float): Current time for external state variable evaluation.
                Only needed if material depends on time (default: 0.0)
              fourierMode (int): Fourier mode (default: -1)
              varc_curr (FieldOnCellsReal): external state variables at current time (default: None)
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
            matr_elem = self.getElasticStiffnessMatrix(
                time, fourierMode, varc_curr, groupOfCells, with_dual
            )

            if assembly:
                matr_asse = AssemblyMatrixDisplacementReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isThermal():
            matr_elem = self.getLinearConductivityMatrix(
                time, fourierMode, varc_curr, groupOfCells, with_dual
            )

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
    def getMassMatrix(self, time=0.0, varc_curr=None, groupOfCells=[], assembly=False):
        """Return the elementary matrices formass matrix depending of the physic.
        Option MASS_MECA or MASS_THER or MASS_ACOU.

        Arguments:
              time (float): Current time for external state variable evaluation.
                Only needed if material depends on time (default: 0.0)
              varc_curr (FieldOnCellsReal): external state variables at current time (default: None)
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it is empty, the full model is used
              assembly (bool): assemble elementary matrix (default: False)
        Returns:
              ElementaryMatrix: elementary mass matrix
        """

        model = self.getPhysicalProblem().getModel()

        if model.isMechanical():
            matr_elem = self.getMechanicalMassMatrix(False, varc_curr, groupOfCells)

            if assembly:
                matr_asse = AssemblyMatrixDisplacementReal(self.getPhysicalProblem())
                matr_asse.addElementaryMatrix(matr_elem)
                matr_asse.assemble()
                return matr_asse

        elif model.isThermal():
            matr_elem = self.getLinearCapacityMatrix(time, varc_curr, groupOfCells)

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

    @profile
    def getVolumetricForces(
        self, time_curr=0.0, time_step=0.0, theta=1, mode=0, varc_curr=None, assembly=True
    ):
        """Return the volumetric forces field

        Arguments:
                time_curr (float): Current time
                time_step (float): Time increment
                theta (float): Theta parameter for time-integration
                mode (int) : fourier mode
                varc_curr (FieldOnCellsReal): external state variables at current time (default: None)
                assembly (bool): assemble if True

        Returns:
                ElementaryVector: elementary volumetric forces vector if assembly=False
                FieldOnNodes: volumetric forces field if assembly=True
        """

        model = self.getPhysicalProblem().getModel()

        if model.isThermal():
            return self.getThermalVolumetricForces(time_curr, time_step, theta, varc_curr, assembly)
        elif model.isMechanical():
            return self.getMechanicalVolumetricForces(
                time_curr, time_step, theta, mode, varc_curr, assembly
            )
        elif model.isAcoustic():
            return self.getAcousticVolumetricForces(assembly)
        else:
            raise RuntimeError("Not implemented")

    @profile
    def getNeumannForces(
        self, time_curr=0.0, time_step=0.0, theta=1, mode=0, varc_curr=None, assembly=True
    ):
        """Return the Neumann forces field

        Arguments:
                time_curr (float): Current time
                time_step (float): Time increment
                theta (float): Theta parameter for time-integration
                mode (int) : fourier mode
                varc_curr (FieldOnCellsReal): external state variables at current time (default: None)
                assembly (bool): assemble if True

        Returns:
                ElementaryVector: elementary Neumann forces vector if assembly=False
                FieldOnNodes: Neumann forces field if assembly=True
        """

        model = self.getPhysicalProblem().getModel()

        if model.isThermal():
            return self.getThermalNeumannForces(time_curr, time_step, theta, assembly)
        elif model.isMechanical():
            return self.getMechanicalNeumannForces(
                time_curr, time_step, theta, mode, varc_curr, assembly
            )
        elif model.isAcoustic():
            return self.getAcousticNeumannForces(assembly)
        else:
            raise RuntimeError("Not implemented")

    @profile
    def getNonLinearNeumannForces(
        self, primal_prev, primal_step, time_prev, time_step, theta, assembly=True
    ):
        """Return the nonlinear Neumann forces field

        Arguments:
                primal_prev : primal solution at the beginning of the time step
                primal_step : incremental primal solution
                time_prev (float): Previous time at the beginning of the time step
                time_step (float): Time increment
                theta (float): Theta parameter for time-integration
                assembly (bool): assemble if True

        Returns:
                ElementaryVector: elementary Neumann forces vector if assembly=False
                FieldOnNodes: Neumann forces field if assembly=True
        """

        model = self.getPhysicalProblem().getModel()

        if model.isThermal():
            return self.getThermalNonLinearNeumannForces(
                primal_prev, primal_step, time_prev, time_step, theta, assembly
            )
        elif model.isMechanical():
            raise RuntimeError("Not implemented")
        elif model.isAcoustic():
            raise RuntimeError("Not implemented")
        else:
            raise RuntimeError("Not implemented")

    @profile
    def getImposedDualBC(self, time_curr=0.0, assembly=True):
        """Return imposed nodal BC field

        Arguments:
                time_curr (float): Current time
                assembly (bool): assemble if True

        Returns:
                ElementaryVector: elementary imposed nodal BC vector if assembly=False
                FieldOnNodes: imposed nodal BC field if assembly=True
        """

        model = self.getPhysicalProblem().getModel()

        if model.isThermal():
            return self.getThermalImposedDualBC(time_curr, assembly)
        elif model.isMechanical():
            return self.getMechanicalImposedDualBC(time_curr, assembly)
        elif model.isAcoustic():
            return self.getAcousticImposedDualBC(assembly)
        else:
            raise RuntimeError("Not implemented")

    @profile
    def getInternalResidual(self, phys_state, scaling=1.0):
        """Compute internal residual R_int(u, Lagr).

            R_int(u, Lagr) = [B^t.Sig(u) + B^t.Lagr, B^t.primal-primal_impo]

        Arguments:
            phys_state (PhysicalState): physical state
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0)

        Returns:
            Residuals: residuals container
            FieldOnCellsReal: internal state variables (VARI_ELGA)
            FieldOnCellsReal: Cauchy stress tensor (SIEF_ELGA)
        """

        # Compute internal forces (B^t.stress)
        _, codret, internVar, stress, r_stress = self.getInternalForces(
            phys_state.primal,
            phys_state.primal_step,
            phys_state.stress,
            phys_state.internVar,
            phys_state.time,
            phys_state.time_step,
            phys_state.getState(-1).externVar,
            phys_state.externVar,
        )

        resi = Residuals()

        resi.resi_stress = r_stress
        r_int = r_stress

        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        phys_pb = self.getPhysicalProblem()
        if phys_pb.getDOFNumbering().useLagrangeDOF():
            primal_curr = phys_state.primal + phys_state.primal_step

            # Compute kinematic forces (B^t.Lagr_curr)
            dualizedBC_forces = self.getDualForces(primal_curr)
            resi.resi_dual = dualizedBC_forces

            # Compute dualized BC (B^t.primal_curr - primal_impo)
            # Compute dualized BC (B^t.primal_curr)
            dualizedBC_disp = self.getDualDisplacement(primal_curr, scaling)

            # Imposed dualized BC (primal_impo)
            time_curr = phys_state.time + phys_state.time_step
            dualizedBC_impo = self.getImposedDualBC(time_curr)

            r_int += dualizedBC_forces + dualizedBC_disp - dualizedBC_impo
        else:
            resi.resi_dual = phys_state.createPrimal(phys_pb, 0.0)

        resi.resi_int = r_int

        return resi, internVar, stress

    @profile
    def getExternalResidual(self, phys_state):
        """Compute external residual R_ext(u, Lagr)

        R_ext(u, Lagr) = [(f(u),v), 0]

        Arguments:
            phys_state (PhysicalState) : physical state

        Returns
            FieldOnNodesReal: external residual field
        """

        # Compute Neumann forces
        neumann_forces = self.getNeumannForces(
            phys_state.time + phys_state.time_step, varc_curr=phys_state.externVar
        )

        volum_forces = self.getVolumetricForces(
            phys_state.time + phys_state.time_step, varc_curr=phys_state.externVar
        )

        resi_ext = neumann_forces + volum_forces

        return resi_ext

    @profile
    def getContactResidual(self, phys_state, contact_manager):
        """Compute contact residual R_cont(u, Lagr)

        Arguments:
            phys_state (PhysicalState) : physical state
            contact_manager (ContactManager) : contact manager

        Returns
            FieldOnNodesReal: contact residual field

        """

        if contact_manager:
            # Compute contact forces
            contact_forces = self.getContactForces(
                contact_manager.getPairingCoordinates(),
                phys_state.primal,
                phys_state.primal_step,
                phys_state.time,
                phys_state.time_step,
                contact_manager.data(),
                contact_manager.coef_cont,
                contact_manager.coef_frot,
            )
        else:
            contact_forces = phys_state.createPrimal(self.getPhysicalProblem(), 0.0)

        return contact_forces

    @profile
    def getResidual(self, phys_state, contact_manager=None, scaling=1.0):
        """Compute R(u, Lagr) = - (Rint(u, Lagr) + Rcont(u, Lagr) - Rext(u, Lagr)).

        This is not the true residual but the opposite.

        Arguments:
            phys_state (PhysicalState) : physical state
            contact_manager (ContactManager) : contact manager
            scaling (float): Scaling factor for Lagrange multipliers (default: 1.0).

        Returns:
            tuple(Residuals, FieldOnCellsReal, FieldOnCellsReal):
            Tuple with residuals, internal state variables (VARI_ELGA),
            Cauchy stress tensor (SIEF_ELGA).
        """

        resi, internVar, stress = self.getInternalResidual(phys_state, scaling)
        resi.resi_ext = self.getExternalResidual(phys_state)
        resi.resi_cont = self.getContactResidual(phys_state, contact_manager)

        # Compute residual
        resi.resi = -(resi.resi_int + resi.resi_cont - resi.resi_ext)

        return resi, internVar, stress

    @profile
    def getInternalTangentMatrix(self, phys_state, matrix_type="TANGENTE", assemble=False):
        """Compute K(u) = d(Rint(u)) / du

        Arguments:
            phys_state (PhysicalState) : physical state
            matrix_type (str): type of matrix used.
            assemble (bool): assemble or not the matrix

        Returns:
            int : error code flag
            ElementaryMatrixDisplacementReal: rigidity matrix.
            ElementaryMatrixDisplacementReal: dual matrix.
        """
        # Main object for discrete computation

        # Compute rigidity matrix
        if matrix_type in ("PRED_ELASTIQUE", "ELASTIQUE"):
            time_curr = phys_state.time + phys_state.time_step
            matr_elem_rigi = self.getLinearStiffnessMatrix(time=time_curr, with_dual=False)
            codret = 0
        elif matrix_type == "PRED_TANGENTE":
            _, codret, matr_elem_rigi = self.getPredictionTangentStiffnessMatrix(
                phys_state.primal,
                phys_state.primal_step,
                phys_state.stress,
                phys_state.internVar,
                phys_state.time,
                phys_state.time_step,
                phys_state.getState(-1).externVar,
                phys_state.externVar,
            )
        elif matrix_type == "TANGENTE":
            _, codret, matr_elem_rigi = self.getTangentStiffnessMatrix(
                phys_state.primal,
                phys_state.primal_step,
                phys_state.stress,
                phys_state.internVar,
                phys_state.time,
                phys_state.time_step,
                phys_state.getState(-1).externVar,
                phys_state.externVar,
            )
        else:
            raise RuntimeError("Matrix not supported: %s" % (matrix_type))

        # Compute dual matrix
        matr_elem_dual = self.getDualStiffnessMatrix()

        if assemble:
            if codret > 0:
                raise IntegrationError("MECANONLINE10_1")

            # Assemble matrix
            jacobian = AssemblyMatrixDisplacementReal(self.getPhysicalProblem())
            jacobian.addElementaryMatrix(matr_elem_rigi)
            jacobian.addElementaryMatrix(matr_elem_dual)

            jacobian.assemble()

            return jacobian

        return codret, matr_elem_rigi, matr_elem_dual

    @profile
    def getContactTangentMatrix(self, phys_state, contact_manager):
        """Compute K(u) = d(Rcont(u) ) / du

        Arguments:
            phys_state (PhysicalState) : physical state
            contact_manager (ContactManager) : contact manager

        Returns:
           ElementaryMatrixDisplacementReal: Contact matrix.
        """
        if contact_manager:
            matr_elem_cont = self.getContactMatrix(
                contact_manager.getPairingCoordinates(),
                phys_state.primal,
                phys_state.primal_step,
                phys_state.time,
                phys_state.time_step,
                contact_manager.data(),
                contact_manager.coef_cont,
                contact_manager.coef_frot,
            )

            return matr_elem_cont

        return None

    @profile
    def getTangentMatrix(
        self, phys_state, matrix_type="TANGENTE", contact_manager=None, assemble=True
    ):
        """Compute tangent matrix for nonlinear problem.
            K(u) = d(Rint(u) - Rext(u)) / du

        Arguments:
            phys_state (PhysicalState) : physical state
            matrix_type (str): type of matrix used.
            contact_manager (ContactManager) : contact manager
            assemble (bool): assemble or not the matrix

        Returns:
            AssemblyMatrixDisplacementReal: Tangent matrix.
        """

        # Compute elementary matrix
        codret, matr_elem_rigi, matr_elem_dual = self.getInternalTangentMatrix(
            phys_state, matrix_type, False
        )
        if codret > 0:
            raise IntegrationError("MECANONLINE10_1")

        matr_elem_cont = self.getContactTangentMatrix(phys_state, contact_manager)

        if assemble:
            # Assemble matrix
            jacobian = AssemblyMatrixDisplacementReal(self.getPhysicalProblem())
            jacobian.addElementaryMatrix(matr_elem_rigi)
            jacobian.addElementaryMatrix(matr_elem_dual)
            jacobian.addElementaryMatrix(matr_elem_cont)

            jacobian.assemble()

            return jacobian

        return matr_elem_rigi, matr_elem_dual, matr_elem_cont
