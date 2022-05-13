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
        .def( "imposedDisplacement", &DiscreteComputation::imposedDisplacement,
              R"(
      Return the imposed displacement assembled vector

      Arguments:
            float: current time

      Returns:
            FieldOnNodes: imposed displacement
        )",
              py::arg( "time" ) )
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
            time_list (list): Vector of times of length 3 (current time, delta_time, parameter)
            externVarField (fieldOnCellsReal): external state variable at current time

      Returns:
            FieldOnNodes: Neumann load vector
        )",
              py::arg( "time_list" ), py::arg( "externVarField" ) = nullptr )
        .def( "createExternalStateVariablesField",
              &DiscreteComputation::createExternalStateVariablesField, R"(
            Create external state variable field

            Arguments:
                  time (float): current time

            Returns:
                  FieldOnCells: field of external state variables at current time
            )",
              py::arg( "time" ) )
        .def( "createTimeField", &DiscreteComputation::createTimeField, R"(
            Create time field

            Arguments:
                  time (float): current time

            Returns:
                  ConstantFieldOnCells: field of current time

            )",
              py::arg( "time" ) )
        .def( "computeExternalStateVariablesLoad",
              &DiscreteComputation::computeExternalStateVariablesLoad, R"(
            Compute load from external state variables

            Arguments:
                  time (float): current time
                  timeField (ConstantFieldOnCell): field with value of current time
                  externVarField (fieldOnCellsReal): external state variable at current time

            Returns:
                  FieldOnNodes: load from external state variables
            )",
              py::arg( "time" ), py::arg( "timeField" ), py::arg( "externVarField" ) )

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
                  time (float): current time

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
                  time (float): current time
                  disp (FieldOnNodes): displacement field at current time

            Returns:
                  FieldOnNodes: incremental imposed displacement vector
        )",
              py::arg( "time" ), py::arg( "disp" ) )
        .def( "elasticStiffnessMatrix", &DiscreteComputation::elasticStiffnessMatrix, R"(
            Return the elementary matrices for elastic Stiffness matrix

            Arguments:
                  time (float): current time (default: 0.0)
                  fourierMode (int): Fourier mode (default: 0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                  ElementaryMatrix: elementary elastic Stiffness matrix
            )",
              py::arg( "time" ) = 0.0, py::arg( "fourierMode" ) = 0,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "externVarField" ) = nullptr )
        .def( "getPhysicalProblem", &DiscreteComputation::getPhysicalProblem, R"(
            Get physical probelm

            Returns:
                  PhysicalProblem: physical problem
            )" )
        .def( "dualStiffnessMatrix", &DiscreteComputation::dualStiffnessMatrix,
              R"(
      Return elementary matrices for dual BC

      Arguments:
            None

      Returns:
            ElementaryMatrix: elementary matrices
        )" )

        .def( "linearConductivityMatrix", &DiscreteComputation::linearConductivityMatrix,
              R"(
      Return the elementary matices for linear thermal matrix

      Arguments:
            time (float): current time
                  fourierMode (int): Fourier mode (default: 0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  externVarField (fieldOnCellsReal): external state variable at current time
      Returns:
            ElementaryMatrix: elementary linear thermal matrices
        )",
              py::arg( "time" ) = 0., py::arg( "delta_time" ) = 0., py::arg( "fourierMode" ) = 0,
              py::arg( "groupOfCells" ) = VectorString(), py::arg( "externVarField" ) = nullptr )

        .def( "massMatrix", &DiscreteComputation::massMatrix, R"(
            Return the elementary matrices for elastic Stiffness matrix

            Arguments:
                  time (float): current time (default: 0.0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                  ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "time" ) = 0., py::arg( "groupOfCells" ) = VectorString(),
              py::arg( "externVarField" ) = nullptr )

        .def( "linearCapacityMatrix", &DiscreteComputation::linearCapacityMatrix, R"(
            Return the elementary matrices for linear Capacity matrix in thermal computation

            Arguments:
                  time (float): current time (default: 0.0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                  ElementaryMatrix: elementary mass matrix
            )",
              py::arg( "time" ) = 0., py::arg( "groupOfCells" ) = VectorString(),
              py::arg( "externVarField" ) = nullptr )

        .def( "dampingMatrix", &DiscreteComputation::dampingMatrix, R"(
            Return the elementary matrices for elastic Stiffness matrix

            Arguments:
                  massMatrix : elementary mass matrix
                  stiffnessMatrix : elementary stiffness matrix
                  time (float): current time (default: 0.0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                      If it empty, the full model is used
                  externVarField (fieldOnCellsReal): external state variable at current time
            Returns:
                  ElementaryMatrix: elementary damping matrix
            )",
              py::arg( "massMatrix" ) = nullptr, py::arg( "stiffnessMatrix" ) = nullptr,
              py::arg( "time" ) = 0., py::arg( "groupOfCells" ) = VectorString(),
              py::arg( "externVarField" ) = nullptr )

        .def( "computeInternalForces", &DiscreteComputation::computeInternalForces,
              R"(
      Compute internal forces (integration of behaviour)

      Arguments:
            displ (FieldOnNodes): displacement field at begin of current time
            displ_incr (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): field of internal state variables at begin of current time
            timeFieldPrev (constantFieldOnCells): time at begin of current time
            timeFieldCurr (constantFieldOnCells): time at end of current time
            groupOfCells (list[str]): compute matrices on given groups of cells.

      Returns:
            tuple (tuple): return code error (FieldOnCells),
            error code flag (integer),
            internal state variables VARI_ELGA (FieldOnCells),
            Cauchy stress SIEF_ELGA (FieldOnCells),
            field of internal forces (FieldOnNodesReal),
        )",
              py::arg( "displ" ), py::arg( "displ_incr" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "timeFieldPrev" ), py::arg( "timeFieldCurr" ),
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "computeTangentStiffnessMatrix", &DiscreteComputation::computeTangentStiffnessMatrix,
              R"(
      Compute jacobian matrix for Newton algorithm

      Arguments:
            displ (FieldOnNodes): displacement field at begin of current time
            displ_incr (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): field of internal state variables at begin of current time
            timeFieldPrev (ConstantFieldOnCells): time at begin of current time
            timeFieldCurr (ConstantFieldOnCells): time at end of current time
            groupOfCells (list[str]): compute matrices on given groups of cells.

      Returns:
            tuple (tuple): return code error (FieldOnCellsLong),
            error code flag (int),
            elementary tangent matrix (ElementaryMatrixDisplacementReal)
        )",
              py::arg( "displ" ), py::arg( "displ_incr" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "timeFieldPrev" ), py::arg( "timeFieldCurr" ),
              py::arg( "groupOfCells" ) = VectorString() )

        .def( "computeTangentPredictionMatrix",
              &DiscreteComputation::computeTangentPredictionMatrix, R"(
      Compute jacobian matrix for Newton algorithm, Euler prediction

      Arguments:
            displ (FieldOnNodes): displacement field at begin of current time
            displ_incr (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): field of internal state variables at begin of current time
            timeFieldPrev (constantFieldOnCells): time at begin of current time
            timeFieldCurr (constantFieldOnCells): time at end of current time
            groupOfCells (list[str]): compute matrices on given groups of cells.

      Returns:
            tuple (tuple): return code error (FieldOnCellsLong),
            error code flag (int),
            elementary tangent matrix (ElementaryMatrixDisplacementReal),
        )",
              py::arg( "displ" ), py::arg( "displ_incr" ), py::arg( "stress" ),
              py::arg( "internVar" ), py::arg( "timeFieldPrev" ), py::arg( "timeFieldCurr" ),
              py::arg( "groupOfCells" ) = VectorString() );
};
