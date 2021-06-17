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

void exportDiscreteProblemToPython() {

    py::class_< DiscreteProblem, DiscreteProblem::DiscreteProblemPtr >( "DiscreteProblem",
                                                                                  py::no_init )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< DiscreteProblem, StudyDescriptionPtr > ) )
        // fake initFactoryPtr: not a DataStructure
        .def( "computeElementaryMechanicalLoadsVector",
              &DiscreteProblem::computeElementaryMechanicalLoadsVector )
        .def( "computeElementaryDirichletVector",
              &DiscreteProblem::computeElementaryDirichletVector )
        .def( "computeElementaryLaplaceVector", &DiscreteProblem::computeElementaryLaplaceVector )
        .def( "computeElementaryNeumannVector", &DiscreteProblem::computeElementaryNeumannVector )
        .def( "computeElementaryStiffnessMatrix",
              &DiscreteProblem::computeElementaryStiffnessMatrix )
        .def( "computeElementaryTangentMatrix", &DiscreteProblem::computeElementaryTangentMatrix )
        .def( "computeElementaryJacobianMatrix",
              &DiscreteProblem::computeElementaryJacobianMatrix )
        .def( "computeDirichletBC", &DiscreteProblem::computeDirichletBC,
              computeDirichletBC_overloads() )
        .def( "computeDOFNumbering", &DiscreteProblem::computeDOFNumbering )
        .def( "computeMechanicalDampingMatrix",
              &DiscreteProblem::computeMechanicalDampingMatrix )
        .def( "computeMechanicalStiffnessMatrix",
              &DiscreteProblem::computeMechanicalStiffnessMatrix )
        .def( "computeMechanicalMassMatrix", &DiscreteProblem::computeMechanicalMassMatrix )
        .def( "getStudyDescription", &DiscreteProblem::getStudyDescription )

        .def( "createBehaviour",
              static_cast< void ( DiscreteProblem::* )( PyObject *, const std::string &,
                                                             const std::string &, const int ) >(
                  &DiscreteProblem::createBehaviour ),
              R"(
Create constant fields on cells for behaviour (COMPOR, CARCRI and MULCOM)

Arguments:
    COMPORTEMENT (list[dict]): keywords as provided to STAT_NON_LINE/COMPORTEMENT
    SIGM_INIT (str): "OUI" if there is an initial stress field
    IMPLEX (str): "OUI" if Implex algorithm is used
    INFO (int): level of verbosity, 1 to have description of behaviour or 0 to be quiet
        )",
              ( py::arg( "self" ), py::arg( "COMPORTEMENT" ), py::arg( "SIGM_INIT" ),
                py::arg( "IMPLEX" ), py::arg( "INFO" ) ) )
        .def( "createBehaviour",
              static_cast< void ( DiscreteProblem::* )( PyObject * ) >(
                  &DiscreteProblem::createBehaviour ),
              R"(
Create constant fields on cells for behaviour (COMPOR, CARCRI and MULCOM)

Arguments:
    COMPORTEMENT (list[dict]): keywords as provided to STAT_NON_LINE/COMPORTEMENT
        )",
              ( py::arg( "self" ), py::arg( "COMPORTEMENT" ) ) );
};
