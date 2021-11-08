/**
 * @file DiscreteComputationInterface.cxx
 * @brief Interface python de DiscreteComputation
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( computeDualizedDirichlet_overloads,
    dualDisplacement, 1, 2 )

void exportDiscreteComputationToPython() {

    py::class_< DiscreteComputation, DiscreteComputation::DiscreteComputationPtr >
                                                ( "DiscreteComputation", py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< DiscreteComputation, PhysicalProblemPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( "imposedDisplacement",
              &DiscreteComputation::imposedDisplacement )
        .def( "dualReaction",
              &DiscreteComputation::dualReaction )
        .def( "dualDisplacement",
              &DiscreteComputation::dualDisplacement,
              computeDualizedDirichlet_overloads() )
        .def( "Neumann", &DiscreteComputation::Neumann )
        .def( "DirichletBC", &DiscreteComputation::DirichletBC )
        .def( "computeElementaryStiffnessMatrix",
              &DiscreteComputation::computeElementaryStiffnessMatrix )
        .def( "computeElementaryTangentMatrix", 
                              &DiscreteComputation::computeElementaryTangentMatrix )
        .def( "computeElementaryJacobianMatrix",
              &DiscreteComputation::computeElementaryJacobianMatrix )
        .def( "computeMechanicalDampingMatrix",
              &DiscreteComputation::computeMechanicalDampingMatrix )
        .def( "computeMechanicalStiffnessMatrix",
              &DiscreteComputation::computeMechanicalStiffnessMatrix )
        .def( "computeMechanicalMassMatrix", &DiscreteComputation::computeMechanicalMassMatrix )
        .def( "getPhysicalProblem", &DiscreteComputation::getPhysicalProblem );
};
