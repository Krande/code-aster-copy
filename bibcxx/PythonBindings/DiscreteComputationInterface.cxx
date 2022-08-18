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
        .def( "imposedDualBC",
              py::overload_cast< const ASTERDOUBLE, const ASTERDOUBLE, const ASTERDOUBLE >(
                  &DiscreteComputation::imposedDualBC, py::const_ ),
              R"(
      Return the imposed nodal BC assembled vector

      Arguments:
            time_value (float): Current time
            time_delta (float): Time increment
            time_theta (float): Theta parameter for integration

      Returns:
            FieldOnNodes: imposed dual field
        )",
              py::arg( "time_value" ), py::arg( "time_delta" ), py::arg( "time_theta" ) )
        .def( "imposedDualBC",
              py::overload_cast< const ASTERDOUBLE >( &DiscreteComputation::imposedDualBC,
                                                      py::const_ ),
              R"(
      Return the imposed nodal BC assembled vector

      Arguments:
            time_value (float): Current time

      Returns:
            FieldOnNodes: imposed dual field
        )",
              py::arg( "time_value" ) )
        .def( "dualReaction", &DiscreteComputation::dualReaction,
              R"(
      Return the imposed displacement assembled vector

      Arguments:
            disp_curr (FieldOnNodes): current displacement vector

      Returns:
            FieldOnNodes: dual reaction vector (B^T*lambda)
        )",
              py::arg( "disp_curr" ) )
        .def( "dualDisplacement", &DiscreteComputation::dualDisplacement,
              R"(
      Return the Dirichlet load vector

      Arguments:
            disp_curr (FieldOnNodes): current displacement vector

      Returns:
            FieldOnNodes: Dirichlet load vector
              )",
              py::arg( "disp_curr" ), py::arg( "scaling" ) = 1.0 )
        .def( "neumann", &DiscreteComputation::neumann,
              R"(
      Return the Neumann load vector

      Arguments:
            time_value (float): Current time
            time_delta (float): Time increment
            time_theta (float): Theta parameter for integration
            externVarField (fieldOnCellsReal): external state variable at current time
            previousPrimalField (fieldOnNodesReal): solution field at previous time

      Returns:
            FieldOnNodes: Neumann load vector
        )",
              py::arg( "time_value" ), py::arg( "time_delta" ), py::arg( "time_theta" ),
              py::arg( "externVarField" ) = nullptr, py::arg( "previousPrimalField" ) = nullptr )
        .def( "createExternalStateVariablesField",
              &DiscreteComputation::createExternalStateVariablesField, R"(
            Create external state variable field

            Arguments:
                  time_value (float): Current time

            Returns:
                  FieldOnCells: field of external state variables at current time
            )",
              py::arg( "time_value" ) )
        .def( "computeExternalStateVariablesLoad",
              &DiscreteComputation::computeExternalStateVariablesLoad, R"(
            Compute load from external state variables

            Arguments:
                  time_value (float): Current time
                  externVarField (fieldOnCellsReal): external state variable at current time

            Returns:
                  FieldOnNodes: load from external state variables
            )",
              py::arg( "time_value" ), py::arg( "externVarField" ) )
        .def( "transientThermalLoad", &DiscreteComputation::transientThermalLoad, R"(
            Compute Transient Thermal Load

            Arguments:
                  time_value (float): Current time
                  time_delta (float): Time increment
                  time_theta (float): Theta parameter for integration
                  externVarField (fieldOnCellsReal): external state variable at current time
                  previousPrimalField (fieldOnNodesReal): solution field at previous time

            Returns:
                  FieldOnNodes: load from external state variables
            )",
              py::arg( "time_value" ), py::arg( "time_delta" ), py::arg( "time_theta" ),
              py::arg( "externVarField" ), py::arg( "previousPrimalField" ) = nullptr )
        .def( "computeExternalStateVariablesReference",
              &DiscreteComputation::computeExternalStateVariablesReference, R"(
            Compute field for external state variables reference value

            Returns:
                  FieldOnCells: field for external state variables reference values
            )" )
        .def( "dirichletBC", &DiscreteComputation::dirichletBC,
              R"(
            Return the imposed displacement vector used to remove imposed DDL

            Arguments:
                  time_value (float): Current time

            Returns:
                  FieldOnNodes: imposed displacement vector
        )",
              py::arg( "time" ) )
        .def( "incrementalDirichletBC", &DiscreteComputation::incrementalDirichletBC,
              R"(
            Return the incremental imposed displacement vector used to remove imposed DDL
            for incremental resolution.

            incr_disp = dirichletBC(time) - disp, with 0.0 for DDL not imposed

            Arguments:
                  time_value (float): Current time
                  disp (FieldOnNodes): displacement field at current time

            Returns:
                  FieldOnNodes: incremental imposed displacement vector
        )",
              py::arg( "time_value" ), py::arg( "disp" ) )
        .def( "elasticStiffnessMatrix", &DiscreteComputation::elasticStiffnessMatrix, R"(
            Return the elementary matrices for elastic Stiffness matrix.
            Option RIGI_MECA.

            Arguments:
                  time_value (float): Current time (default: 0.0)
                  fourierMode (int): Fourier mode (default: -1)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                  ElementaryMatrix: elementary elastic Stiffness matrix
            )",
              py::arg( "time_value" ) = 0.0, py::arg( "fourierMode" ) = -1,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "externVarField" ) = nullptr )
        .def( "getPhysicalProblem", &DiscreteComputation::getPhysicalProblem, R"(
            Get physical probelm

            Returns:
                  PhysicalProblem: physical problem
            )" )
        .def( "dualStiffnessMatrix", &DiscreteComputation::dualStiffnessMatrix,
              R"(
            Return elementary matrices for dual mechanical BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )

        .def( "dualMobilityMatrix", &DiscreteComputation::dualMobilityMatrix,
              R"(
            Return elementary matrices for dual acoustic BC

            Returns:
                ElementaryMatrix: elementary matrices
        )" )

        .def( "linearConductivityMatrix", &DiscreteComputation::linearConductivityMatrix,
              R"(
            Return the elementary matices for linear thermal matrix.
            Option RIGI_THER.

            Arguments:
                time_value (float): Current time
                time_delta (float): Time increment
                time_theta (float): Theta parameter for integration
                fourierMode (int): Fourier mode (default: -1)
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                ElementaryMatrix: elementary linear thermal matrices
        )",
              py::arg( "time_value" ), py::arg( "time_delta" ), py::arg( "time_theta" ),
              py::arg( "fourierMode" ) = 0, py::arg( "groupOfCells" ) = VectorString(),
              py::arg( "externVarField" ) = nullptr )

        .def( "linearMobilityMatrix", &DiscreteComputation::linearMobilityMatrix,
              R"(
            Return the elementary matices for linear mobility acoustic matrix
            Option RIGI_ACOU.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
            Returns:
                ElementaryMatrix: elementary linear acoustic matrices
        )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "massMatrix", &DiscreteComputation::massMatrix, R"(
            Return the elementary matrices for mechanical mass matrix
            Option MASS_MECA.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "externVarField" ) = nullptr )

        .def( "compressibilityMatrix", &DiscreteComputation::compressibilityMatrix, R"(
            Return the elementary matrices for compressibility acoustic matrix.
            Option MASS_ACOU.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "linearCapacityMatrix", &DiscreteComputation::linearCapacityMatrix, R"(
            Return the elementary matrices for linear Capacity matrix in thermal computation.
            Option MASS_THER.

            Arguments:
                time_value (float): Current time
                time_delta (float): Time increment
                time_theta (float): Theta parameter for integration
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "time_value" ), py::arg( "time_delta" ), py::arg( "time_theta" ),
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "externVarField" ) = nullptr )

        .def( "dampingMatrix", &DiscreteComputation::dampingMatrix, R"(
            Return the elementary matrices for damping matrix.
            Option AMOR_MECA.

            Arguments:
                massMatrix : elementary mass matrix
                stiffnessMatrix : elementary stiffness matrix
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                ElementaryMatrixReal: elementary damping matrix
            )",
              py::arg( "massMatrix" ) = nullptr, py::arg( "stiffnessMatrix" ) = nullptr,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "externVarField" ) = nullptr )

        .def( "impedanceMatrix", &DiscreteComputation::impedanceMatrix, R"(
            Return the elementary matrices for impedance (acoustic) damping matrix.
            Option AMOR_ACOU.

            Returns:
                ElementaryMatrixReal: elementary damping matrix
            )" )

        .def( "complexStiffnessMatrix", &DiscreteComputation::complexStiffnessMatrix, R"(
            Return the elementary matrices for viscoelastic Stiffness matrix.
            Option RIGI_MECA_HYST.

            Arguments:
                stiffnessMatrix : elementary stiffness matrix
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                ElementaryMatrixComplex: elementary viscoelastic rigidity matrix
            )",
              py::arg( "stiffnessMatrix" ), py::arg( "groupOfCells" ) = VectorString(),
              py::arg( "externVarField" ) = nullptr )

        .def( "gyroscopicStiffnessMatrix", &DiscreteComputation::gyroscopicStiffnessMatrix, R"(
            Return the elementary matrices for gyroscopic Stiffness matrix.
            Option RIGI_GYRO.

            Arguments:
                groupOfCells (list[str]): compute matrices on given groups of cells.
                    If it empty, the full model is used
                externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                ElementaryMatrixReal: elementary gyroscopic rigidity matrix
            )",
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "externVarField" ) = nullptr )

        .def( "computeInternalForces", &DiscreteComputation::computeInternalForces,
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

        .def( "computeTangentStiffnessMatrix", &DiscreteComputation::computeTangentStiffnessMatrix,
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
        .def( "contactForces", &DiscreteComputation::contactForces, R"(
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
        .def( "contactMatrix", &DiscreteComputation::contactMatrix, R"(
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
