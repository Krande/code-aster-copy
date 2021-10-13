/**
 * @file DiscreteProblemInterface.cxx
 * @brief Interface python de DiscreteProblem
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
#include "PythonBindings/DiscreteProblemInterface.h"
#include <PythonBindings/factory.h>

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( computeDirichletBC_overloads, computeDirichletBC, 2, 2 )
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( computeElementaryDualizedDirichletVector_overloads,
    computeElementaryDualizedDirichletVector, 1, 2 )
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( computeDualizedDirichlet_overloads,
    computeDualizedDirichlet, 2, 3 )

void exportDiscreteProblemToPython() {

    py::class_< DiscreteProblem, DiscreteProblem::DiscreteProblemPtr >( "DiscreteProblem",
                                                                                  py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< DiscreteProblem, PhysicalProblemPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( "computeElementaryMechanicalLoadsVector",
              &DiscreteProblem::computeElementaryMechanicalLoadsVector )
        .def( "computeElementaryDirichletVector",
              &DiscreteProblem::computeElementaryDirichletVector )
        .def( "computeDirichlet",
              &DiscreteProblem::computeDirichlet )
        .def( "computeElementaryDirichletReactionVector",
              &DiscreteProblem::computeElementaryDirichletReactionVector )
        .def( "computeDirichletReaction",
              &DiscreteProblem::computeDirichletReaction )
        .def( "computeElementaryDualizedDirichletVector",
              &DiscreteProblem::computeElementaryDualizedDirichletVector,
              computeElementaryDualizedDirichletVector_overloads() )
        .def( "computeDualizedDirichlet",
              &DiscreteProblem::computeDualizedDirichlet,
              computeDualizedDirichlet_overloads() )
        .def( "computeElementaryLaplaceVector", &DiscreteProblem::computeElementaryLaplaceVector )
        .def( "computeElementaryNeumannVector", &DiscreteProblem::computeElementaryNeumannVector )
        .def( "computeNeumann", &DiscreteProblem::computeNeumann )
        .def( "computeElementaryStiffnessMatrix",
              &DiscreteProblem::computeElementaryStiffnessMatrix )
        .def( "computeElementaryTangentMatrix", &DiscreteProblem::computeElementaryTangentMatrix )
        .def( "computeElementaryJacobianMatrix",
              &DiscreteProblem::computeElementaryJacobianMatrix )
        .def( "computeDirichletBC", &DiscreteProblem::computeDirichletBC,
              computeDirichletBC_overloads() )
        .def( "computeMechanicalDampingMatrix",
              &DiscreteProblem::computeMechanicalDampingMatrix )
        .def( "computeMechanicalStiffnessMatrix",
              &DiscreteProblem::computeMechanicalStiffnessMatrix )
        .def( "computeMechanicalMassMatrix", &DiscreteProblem::computeMechanicalMassMatrix )
        .def( "getPhysicalProblem", &DiscreteProblem::getPhysicalProblem );
};
