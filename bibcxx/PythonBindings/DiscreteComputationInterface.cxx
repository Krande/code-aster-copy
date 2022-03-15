/**
 * @file DiscreteComputationInterface.cxx
 * @brief Interface python de DiscreteComputation
 * @author Nicolas Sellenet
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

/* person_in_charge: nicolas.sellenet at edf.fr */

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
            time (list): Vector of times of length 3 (current time, delta_time, parameter)

      Returns:
            FieldOnNodes: Neumann load vector
        )",
              py::arg( "time" ) )
        .def( "externalStateVariables", &DiscreteComputation::externalStateVariables,
              R"(
      Return the external State Variables load vector

      Arguments:
            time (float): current time

      Returns:
            FieldOnNodes: external State Variables load vector
        )",
              py::arg( "time" ) )
        .def( "dirichletBC", &DiscreteComputation::dirichletBC,
              R"(
      Return the imposed displacement vector used to remove imposed DDL

      Arguments:
            float: current time

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
            Return the elementary matices for elastic Stiffness matrix

            Arguments:
                  time (float): current time (default: 0.0)
                  fourierMode (int): Fourier mode (default: 0)
                  groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it empty, the full model is used (default: [])

            Returns:
                  ElementaryMatrix: elementary elastic Stiffness matrices
            )",
              py::arg( "time" ) = 0.0, py::arg( "fourierMode" ) = 0,
              py::arg( "groupOfCells" ) = VectorString() )
        .def( "getPhysicalProblem", &DiscreteComputation::getPhysicalProblem )
        .def( "dualStiffnessMatrix", &DiscreteComputation::dualStiffnessMatrix,
              R"(
      Return elementary matrices for dual BC

      Returns:
            ElementaryMatrix: elementary matrices
        )" )
        .def( "massMatrix", &DiscreteComputation::massMatrix, py::arg( "time" ) = 0. );
};
