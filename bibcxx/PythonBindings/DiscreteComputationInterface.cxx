/**
 * @file DiscreteComputationInterface.cxx
 * @brief Interface python de DiscreteComputation
 * @section LICENCE
 *   Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
 *
 *   This file is part of Code_Aster.
 *
 *   Code_Aster is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Code_Aster is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "PythonBindings/DiscreteComputationInterface.h"

#include "aster_pybind.h"

void exportDiscreteComputationToPython( py::module_ &mod ) {

    py::class_< DiscreteComputation, DiscreteComputation::DiscreteComputationPtr >(
        mod, "DiscreteComputation" )
        .def( py::init( &initFactoryPtr< DiscreteComputation, PhysicalProblemPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( "getImposedDualBC",
              py::overload_cast< const ASTERDOUBLE, const ASTERDOUBLE, const ASTERDOUBLE >(
                  &DiscreteComputation::getImposedDualBC, py::const_ ),
              R"(
      Return the imposed nodal BC assembled vector

      Arguments:
            time (float): Current time
            time_step (float): Time increment
            theta (float): Theta parameter for integration

      Returns:
            FieldOnNodes: imposed dual field
        )",
              py::arg( "time" ), py::arg( "time_step" ), py::arg( "theta" ) )
        .def( "getImposedDualBC",
              py::overload_cast< const ASTERDOUBLE >( &DiscreteComputation::getImposedDualBC,
                                                      py::const_ ),
              R"(
      Return the imposed nodal BC assembled vector

      Arguments:
            time (float): Current time

      Returns:
            FieldOnNodes: imposed dual field
        )",
              py::arg( "time" ) )
        .def( "getDualForces", &DiscreteComputation::getDualForces,
              R"(
      Return the imposed displacement assembled vector

      Arguments:
            disp_curr (FieldOnNodes): current displacement vector

      Returns:
            FieldOnNodes: dual reaction vector (B^T*lambda)
        )",
              py::arg( "disp_curr" ) )
        .def( "getDualDisplacement", &DiscreteComputation::getDualDisplacement,
              R"(
      Return the Dirichlet load vector

      Arguments:
            disp_curr (FieldOnNodes): current displacement vector

      Returns:
            FieldOnNodes: Dirichlet load vector
              )",
              py::arg( "disp_curr" ), py::arg( "scaling" ) = 1.0 )
        .def( "getNeumannForces", &DiscreteComputation::getNeumannForces,
              R"(
      Return the Neumann forces vector

      Arguments:
            time (float): Current time
            time_step (float): Time increment
            theta (float): Theta parameter for time-integration
            previousPrimalField (fieldOnNodesReal): solution field at previous time

      Returns:
            FieldOnNodes: Neumann forces vector
        )",
              py::arg( "time" ), py::arg( "time_step" ), py::arg( "theta" ),
              py::arg( "previousPrimalField" ) = nullptr )
        .def( "getExternalStateVariablesForces",
              &DiscreteComputation::getExternalStateVariablesForces, R"(
            Compute load from external state variables

            Arguments:
                  time (float): Current time

            Returns:
                  FieldOnNodes: load from external state variables
            )",
              py::arg( "time" ) )
        .def( "getTransientThermalForces", &DiscreteComputation::getTransientThermalForces, R"(
            Compute Transient Thermal Load

            Arguments:
                  time (float): Current time
                  time_step (float): Time increment
                  theta (float): Theta parameter for integration
                  previousPrimalField (fieldOnNodesReal): solution field at previous time

            Returns:
                  FieldOnNodes: load from external state variables
            )",
              py::arg( "time" ), py::arg( "time_step" ), py::arg( "theta" ),
              py::arg( "previousPrimalField" ) = nullptr )

        .def( "getDirichletBC", &DiscreteComputation::getDirichletBC,
              R"(
            Return the imposed displacement vector used to remove imposed DDL

            Arguments:
                  time (float): Current time

            Returns:
                  FieldOnNodes: imposed displacement vector
        )",
              py::arg( "time" ) )
        .def( "getIncrementalDirichletBC", &DiscreteComputation::getIncrementalDirichletBC,
              R"(
            Return the incremental imposed displacement vector used to remove imposed DDL
            for incremental resolution.

            incr_disp = getDirichletBC(time) - disp, with 0.0 for DDL not imposed

            Arguments:
                  time (float): Current time
                  disp (FieldOnNodes): displacement field at current time

            Returns:
                  FieldOnNodes: incremental imposed displacement vector
        )",
              py::arg( "time" ), py::arg( "disp" ) )
        .def( "getElasticStiffnessMatrix", &DiscreteComputation::getElasticStiffnessMatrix, R"(
            Return the elementary matrices for elastic Stiffness matrix.
            Option RIGI_MECA.

            Arguments:
                  time (float): Current time for external state variavle evaluation (default: 0.0)
                  fourierMode (int): Fourier mode (default: -1)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  with_dual (bool): compute dual terms or not (default: True)
            Returns:
                  ElementaryMatrix: elementary elastic Stiffness matrix
            )",
              py::arg( "time" ) = 0.0, py::arg( "fourierMode" ) = -1,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "with_dual" ) = true )

        .def( "getFluidStructureStiffnessMatrix",
              &DiscreteComputation::getFluidStructureStiffnessMatrix,
              R"(
            Return the elementary matrices for fluid-structure stiffness matrix.
            Option RIGI_FLUI_STRUC.

            Arguments:
                  time (float): Current time for external state variavle evaluation (default: 0.0)
                  fourierMode (int): Fourier mode (default: -1)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
            Returns:
                  ElementaryMatrixReal: elementary fluid-structure Stiffness matrix
            )",
              py::arg( "time" ) = 0.0, py::arg( "fourierMode" ) = -1,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getFluidStructureMassMatrix", &DiscreteComputation::getFluidStructureMassMatrix,
              R"(
            Return the elementary matrices for fluid-structure mass matrix.
            Option MASS_FLUI_STRUC.

            Arguments:
                  time (float): Current time for external state variavle evaluation (default: 0.0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
            Returns:
                  ElementaryMatrixReal: elementary fluid-structure mass matrix
            )",
              py::arg( "time" ) = 0.0, py::arg( "groupOfCells" ) = VectorString() )

        .def( "getPhysicalProblem", &DiscreteComputation::getPhysicalProblem, R"(
            Get physical probelm

            Returns:
                  PhysicalProblem: physical problem
            )" )
        .def( "getDualElasticStiffnessMatrix", &DiscreteComputation::getDualElasticStiffnessMatrix,
              R"(
            Return elementary matrices for dual mechanical BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )

        .def( "getDualLinearMobilityMatrix", &DiscreteComputation::getDualLinearMobilityMatrix,
              R"(
            Return elementary matrices for dual acoustic BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )
        .def( "getDualLinearConductivityMatrix",
              &DiscreteComputation::getDualLinearConductivityMatrix,
              R"(
            Return elementary matrices for dual thermal BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )
        .def( "getLinearConductivityMatrix", &DiscreteComputation::getLinearConductivityMatrix,
              R"(
            Return the elementary matices for linear thermal matrix.
            Option RIGI_THER.

            Arguments:
                  time (float): Current time
                  fourierMode (int): Fourier mode (default: -1)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                  with_dual (bool): compute dual terms or not (default: True)
            Returns:
                  ElementaryMatrix: elementary linear thermal matrices
        )",
              py::arg( "time" ), py::arg( "fourierMode" ) = 0,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "with_dual" ) = true )

        .def( "getExchangeThermalMatrix", &DiscreteComputation::getExchangeThermalMatrix,
              R"(
            Return the elementary matices for exhange thermal matrix.

            Arguments:
                time (float): Current time
            Returns:
                ElementaryMatrix: elementary exchange thermal matrices
        )",
              py::arg( "time" ) )

        .def( "getLinearMobilityMatrix", &DiscreteComputation::getLinearMobilityMatrix,
              R"(
            Return the elementary matices for linear mobility acoustic matrix
            Option RIGI_ACOU.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                with_dual (bool): compute dual terms or not (default: True)

            Returns:
                ElementaryMatrix: elementary linear acoustic matrices
        )",
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "with_dual" ) = true )

        .def( "getMechanicalMassMatrix", &DiscreteComputation::getMechanicalMassMatrix, R"(
            Return the elementary matrices for mechanical mass matrix
            Option MASS_MECA.

            Arguments:
                diagonal (bool) : True for diagonal mass matrix else False.
                time (float): Current time for external state variavle evaluation (default: 0.0)
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "diagonal" ), py::arg( "time" ) = 0.0,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getCompressibilityMatrix", &DiscreteComputation::getCompressibilityMatrix, R"(
            Return the elementary matrices for compressibility acoustic matrix.
            Option MASS_ACOU.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getLinearCapacityMatrix", &DiscreteComputation::getLinearCapacityMatrix, R"(
            Return the elementary matrices for linear Capacity matrix in thermal computation.
            Option MASS_THER.

            Arguments:
                time (float): current time to evaluate rho_cp
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "time" ), py::arg( "groupOfCells" ) = VectorString() )

        .def( "getMechanicalDampingMatrix", &DiscreteComputation::getMechanicalDampingMatrix, R"(
            Return the elementary matrices for damping matrix.
            Option AMOR_MECA.

            Arguments:
                getMechanicalMassMatrix : elementary mass matrix
                stiffnessMatrix : elementary stiffness matrix
                time (float): Current time for external state variavle evaluation (default: 0.0)
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixReal: elementary damping matrix
            )",
              py::arg( "getMechanicalMassMatrix" ) = nullptr,
              py::arg( "stiffnessMatrix" ) = nullptr, py::arg( "time" ) = 0.0,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getImpedanceMatrix", &DiscreteComputation::getImpedanceMatrix, R"(
            Return the elementary matrices for impedance (acoustic) damping matrix.
            Option AMOR_ACOU.

            Returns:
                ElementaryMatrixReal: elementary damping matrix
            )" )

        .def( "getImpedanceBoundaryMatrix", &DiscreteComputation::getImpedanceBoundaryMatrix, R"(
            Return the elementary matrices for impedance (mechanical) matrix.
            Option IMPE_MECA.

            Returns:
                ElementaryMatrixReal: impedance mechanical matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getImpedanceWaveMatrix", &DiscreteComputation::getImpedanceWaveMatrix, R"(
            Return the elementary matrices for impedance (mechanical) matrix
            from an harmonic wave.
            Option ONDE_FLUI.

            Returns:
                ElementaryMatrixReal: impedance wave matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getHystereticStiffnessMatrix", &DiscreteComputation::getHystereticStiffnessMatrix,
              R"(
            Return the elementary matrices for viscoelastic Stiffness matrix.
            Option RIGI_MECA_HYST.

            Arguments:
                stiffnessMatrix : elementary stiffness matrix
                time (float): Current time for external state variavle evaluation (default: 0.0)
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixComplex: elementary viscoelastic rigidity matrix
            )",
              py::arg( "stiffnessMatrix" ), py::arg( "time" ) = 0.0,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getGeometricStiffnessMatrix", &DiscreteComputation::getGeometricStiffnessMatrix, R"(
            Return the elementary matrices for geometric Stiffness matrix.
            Option RIGI_MECA_HYST.

            Arguments:
                sief_elga (FieldOnCellsReal) : stress at Gauss points
                strx_elga (FieldOnCellsReal) : stress at Gauss points for structural element
                displ (FieldOnNodesReal) : displacement field
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixComplex: elementary geometric rigidity matrix
            )",
              py::arg( "sief_elga" ), py::arg( "strx_elga" ) = nullptr,
              py::arg( "displ" ) = nullptr, py::arg( "modeFourier" ) = -1,
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getRotationalStiffnessMatrix", &DiscreteComputation::getRotationalStiffnessMatrix,
              R"(
            Return the elementary matrices for rotational Stiffness matrix.
            Option RIGI_ROTA.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixReal: elementary rotational rigidity matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getGyroscopicStiffnessMatrix", &DiscreteComputation::getGyroscopicStiffnessMatrix,
              R"(
            Return the elementary matrices for gyroscopic Stiffness matrix.
            Option RIGI_GYRO.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixReal: elementary gyroscopic rigidity matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getGyroscopicDampingMatrix", &DiscreteComputation::getGyroscopicDampingMatrix, R"(
            Return the elementary matrices for gyroscopic damping matrix.
            Option MECA_GYRO.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrixReal: elementary gyroscopic damping matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getInternalForces", &DiscreteComputation::getInternalForces,
              R"(
            Compute internal forces (integration of behaviour)

            Arguments:
                displ (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                stress (FieldOnCells): field of stress at begin of current time
                internVar (FieldOnCells): field of internal state variables at begin of current time
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                groupOfCells (list[str]): compute matrices on given groups of cells.

            Returns:
                tuple (tuple): return code error (FieldOnCells),
                error code flag (integer),
                internal state variables VARI_ELGA (FieldOnCells),
                Cauchy stress SIEF_ELGA (FieldOnCells),
                field of internal forces (FieldOnNodesReal),
        )",
              py::arg( "displ" ), py::arg( "displ_step" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "time_prev" ), py::arg( "time_step" ),
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "getTangentStiffnessMatrix", &DiscreteComputation::getTangentStiffnessMatrix,
              R"(
            Compute jacobian matrix for Newton algorithm

            Arguments:
                displ (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                stress (FieldOnCells): field of stress at begin of current time
                internVar (FieldOnCells): field of internal state variables at begin of current time
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                groupOfCells (list[str]): compute matrices on given groups of cells.

            Returns:
                tuple (tuple): return code error (FieldOnCellsLong),
                error code flag (int),
                elementary tangent matrix (ElementaryMatrixDisplacementReal)
        )",
              py::arg( "displ" ), py::arg( "displ_step" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "time_prev" ), py::arg( "time_step" ),
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "computeTangentPredictionMatrix",
              &DiscreteComputation::computeTangentPredictionMatrix, R"(
            Compute jacobian matrix for Newton algorithm, Euler prediction

            Arguments:
                displ (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                stress (FieldOnCells): field of stress at begin of current time
                internVar (FieldOnCells): field of internal state variables at begin of current time
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                groupOfCells (list[str]): compute matrices on given groups of cells.

            Returns:
                tuple (tuple): return code error (FieldOnCellsLong),
                error code flag (int),
                elementary tangent matrix (ElementaryMatrixDisplacementReal),
        )",
              py::arg( "displ" ), py::arg( "displ_step" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "time_prev" ), py::arg( "time_step" ),
              py::arg( "groupOfCells" ) = VectorString() )
        .def( "getContactForces", &DiscreteComputation::getContactForces, R"(
            Compute contact and friction forces

            Arguments:
                geom (MeshCoordinatesField): coordinates of mesh used to compute normal
                displ (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                data (FieldOnCellsReal): contact data
                coef_cont (FeildOnNodesReal) : contact coefficient
                coef_frot (FeildOnNodesReal) : friction coefficient

            Returns:
                FieldOnNodesReal: contact and friction forces

        )",
              py::arg( "geom" ), py::arg( "displ" ), py::arg( "displ_step" ),
              py::arg( "time_prev" ), py::arg( "time_step" ), py::arg( "data" ),
              py::arg( "coef_cont" ), py::arg( "coef_frot" ) )
        .def( "getContactMatrix", &DiscreteComputation::getContactMatrix, R"(
            Compute contact matrix

            Arguments:
                geom (MeshCoordinatesField): coordinates of mesh used to compute normal
                displ (FieldOnNodes): displacement field at begin of current time
                displ_step (FieldOnNodes): field of increment of displacement
                time_prev (float): time at begin of the step
                time_curr (float): delta time between begin and end of the step
                data (FieldOnCellsReal): contact data
                coef_cont (FeildOnNodesReal) : contact coefficient
                coef_frot (FeildOnNodesReal) : friction coefficient

            Returns:
                ElementaryMatrixDisplacementReal: contact and friction elementary matrix
        )",
              py::arg( "geom" ), py::arg( "displ" ), py::arg( "displ_step" ),
              py::arg( "time_prev" ), py::arg( "time_step" ), py::arg( "data" ),
              py::arg( "coef_cont" ), py::arg( "coef_frot" ) );
};
