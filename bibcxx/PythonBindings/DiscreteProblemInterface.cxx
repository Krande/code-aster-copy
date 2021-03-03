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
#include <PythonBindings/factory.h>
#include "PythonBindings/DiscreteProblemInterface.h"

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( buildDirichletBC_overloads, buildDirichletBC, 2, 2 )

void exportDiscreteProblemToPython() {

    py::class_< DiscreteProblemClass, DiscreteProblemClass::DiscreteProblemPtr >(
        "DiscreteProblem", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< DiscreteProblemClass, StudyDescriptionPtr >))
        // fake initFactoryPtr: not a DataStructure
        .def( "buildElementaryMechanicalLoadsVector",
              &DiscreteProblemClass::buildElementaryMechanicalLoadsVector )
        .def( "buildElementaryDirichletVector",
              &DiscreteProblemClass::buildElementaryDirichletVector )
        .def( "buildElementaryLaplaceVector",
              &DiscreteProblemClass::buildElementaryLaplaceVector )
        .def( "buildElementaryNeumannVector",
              &DiscreteProblemClass::buildElementaryNeumannVector )
        .def( "buildElementaryStiffnessMatrix",
              &DiscreteProblemClass::buildElementaryStiffnessMatrix )
        .def( "buildElementaryTangentMatrix",
              &DiscreteProblemClass::buildElementaryTangentMatrix )
        .def( "buildElementaryJacobianMatrix",
              &DiscreteProblemClass::buildElementaryJacobianMatrix )
        .def( "buildDirichletBC", &DiscreteProblemClass::buildDirichletBC,
              buildDirichletBC_overloads() )
        .def( "computeDOFNumbering", &DiscreteProblemClass::computeDOFNumbering )
        .def( "computeMechanicalDampingMatrix",
              &DiscreteProblemClass::computeMechanicalDampingMatrix )
        .def( "computeMechanicalStiffnessMatrix",
              &DiscreteProblemClass::computeMechanicalStiffnessMatrix )
        .def( "computeMechanicalMassMatrix", &DiscreteProblemClass::computeMechanicalMassMatrix )
        .def( "getStudyDescription", &DiscreteProblemClass::getStudyDescription )
        .def( "createBehaviour", &DiscreteProblemClass::createBehaviour,
               R"(
Create maps for behaviour (COMPOR, CARCRI and MULCOM)

Arguments:
    initialState: set 1 if there is an initial stress
    implex: set 1 if using Implex algorithm
    verbosity: level of verbosity, 1 to have description of behaviour

Returns:
    nothing
        )");
};
