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

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/DiscreteComputationInterface.h"

#include <PythonBindings/factory.h>

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( computeDualizedDirichlet_overloads, dualDisplacement, 1, 2 )

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( computeMassMatrix_overloads, massMatrix, 0, 1 )

void exportDiscreteComputationToPython() {

    py::class_< DiscreteComputation, DiscreteComputation::DiscreteComputationPtr >(
        "DiscreteComputation", py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< DiscreteComputation, PhysicalProblemPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( "imposedDisplacement", &DiscreteComputation::imposedDisplacement,
              R"(
      Return the imposed displacement assembled vector

      Arguments:
            float: current time

      Returns:
            FieldOnNodes: imposed displacement
        )",
              ( py::arg( "self" ), py::arg( "time" ) ) )
        .def( "dualReaction", &DiscreteComputation::dualReaction,
              R"(
      Return the imposed displacement assembled vector

      Arguments:
            disp_curr (FieldOnNodes): current displacement vector

      Returns:
            FieldOnNodes: dual reaction vector (B^T*lambda)
        )",
              ( py::arg( "self" ), py::arg( "disp_curr" ) ) )
        .def( "dualDisplacement", &DiscreteComputation::dualDisplacement,
              computeDualizedDirichlet_overloads() )
        .def( "neumann", &DiscreteComputation::neumann,
              R"(
      Return the Neumann load vector

      Arguments:
            time (list): Vector of times of length 3 (current time, delta_time, parameter)

      Returns:
            FieldOnNodes: Neumann load vector
        )",
              ( py::arg( "self" ), py::arg( "time" ) ) )
        .def( "externalStateVariables", &DiscreteComputation::externalStateVariables,
              R"(
      Return the external State Variables load vector

      Arguments:
            time (float): current time

      Returns:
            FieldOnNodes: external State Variables load vector
        )",
              ( py::arg( "self" ), py::arg( "time" ) ) )
        .def( "dirichletBC", &DiscreteComputation::dirichletBC,
              R"(
      Return the imposed displacement vector used to remove imposed DDL

      Arguments:
            float: current time

      Returns:
            FieldOnNodes: imposed displacement vector
        )",
              ( py::arg( "self" ), py::arg( "time" ) ) )
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
              ( py::arg( "self" ), py::arg( "time" ), py::arg( "disp" ) ) )
        .def( "elasticStiffnessMatrix_", &DiscreteComputation::elasticStiffnessMatrix_,
              R"(
      Return the elementary matices for elastic Stiffness matrix

      Arguments:
            time (float): current time
            fourierMode (int): Fourier mode (>=0)
            groupOfCells (list[str]): compute matrices on given groups of cells.
            If it empty, the full model is used

      Returns:
            ElementaryMatrix: elementary elastic Stiffness matrices
        )",
              ( py::arg( "self" ), py::arg( "time" ), py::arg( "fourierMode" ),
                py::arg( "groupOfCells" ) ) )
        .def( "getPhysicalProblem", &DiscreteComputation::getPhysicalProblem )
        .def( "dualStiffnessMatrix", &DiscreteComputation::dualStiffnessMatrix,
              R"(
      Return elementary matrices for dual BC

      Arguments:
            None

      Returns:
            ElementaryMatrix: elementary matrices
        )" )
        .def( "massMatrix", &DiscreteComputation::massMatrix, computeMassMatrix_overloads() );
};
