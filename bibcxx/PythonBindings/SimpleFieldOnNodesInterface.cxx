/**
 * @file SimpleFieldOnNodesInterface.cxx
 * @brief Interface python de SimpleFieldOnNodes
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2020  EDF R&D                www.code-aster.org
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
#include "PythonBindings/DataStructureInterface.h"
#include "PythonBindings/SimpleFieldOnNodesInterface.h"

void exportSimpleFieldOnNodesToPython() {
    py::class_< SimpleFieldOnNodesRealClass, SimpleFieldOnNodesRealPtr,
                py::bases< DataStructure > >( "SimpleFieldOnNodesReal", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< SimpleFieldOnNodesRealClass >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< SimpleFieldOnNodesRealClass, std::string >))
        .def( "getValue", &SimpleFieldOnNodesRealClass::getValue,
              py::return_value_policy< py::return_by_value >())
        .def( "getValues", &SimpleFieldOnNodesRealClass::getValues,
              py::return_value_policy< py::return_by_value >(), R"(
Return a list of the field values as (x1, y1, z1, x2, y2, z2...)

Returns:
    list[float]: List of values of the field on nodes.
        )", ( py::arg("self" )))
        .def( "getNumberOfComponents", &SimpleFieldOnNodesRealClass::getNumberOfComponents )
        .def( "getNumberOfNodes", &SimpleFieldOnNodesRealClass::getNumberOfNodes )
        .def( "getNameOfComponents", &SimpleFieldOnNodesRealClass::getNameOfComponents )
        .def( "getNameOfComponent", &SimpleFieldOnNodesRealClass::getNameOfComponent )
        .def( "updateValuePointers", &SimpleFieldOnNodesRealClass::updateValuePointers );
    py::class_< SimpleFieldOnNodesComplexClass, SimpleFieldOnNodesComplexPtr,
                py::bases< DataStructure > >( "SimpleFieldOnNodesComplex", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< SimpleFieldOnNodesComplexClass >))
        .def( "__init__", py::make_constructor(
                              &initFactoryPtr< SimpleFieldOnNodesComplexClass, std::string >))
        .def( "getValue", &SimpleFieldOnNodesComplexClass::getValue,
              py::return_value_policy< py::return_by_value >() )
        .def( "getValues", &SimpleFieldOnNodesComplexClass::getValues,
              py::return_value_policy< py::return_by_value >() )
        .def( "getNumberOfComponents", &SimpleFieldOnNodesComplexClass::getNumberOfComponents )
        .def( "getNumberOfNodes", &SimpleFieldOnNodesComplexClass::getNumberOfNodes )
        .def( "getNameOfComponents", &SimpleFieldOnNodesComplexClass::getNameOfComponents )
        .def( "getNameOfComponent", &SimpleFieldOnNodesComplexClass::getNameOfComponent )
        .def( "updateValuePointers", &SimpleFieldOnNodesComplexClass::updateValuePointers );
};
