/**
 * @file BehaviourPropertyInterface.cxx
 * @brief Interface python de BehaviourProperty
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

#include <boost/python.hpp>

namespace py = boost::python;
#include <PythonBindings/factory.h>
#include "PythonBindings/BehaviourPropertyInterface.h"

void exportBehaviourPropertyToPython() {

    py::class_< BehaviourProperty, BehaviourProperty::BehaviourPropertyPtr,
            py::bases< DataStructure > >( "BehaviourProperty", py::no_init )
        .def( "__init__", py::make_constructor(&initFactoryPtr< BehaviourProperty >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< BehaviourProperty, ModelPtr, MaterialFieldPtr >))
        .def("getModel", &BehaviourProperty::getModel, R"(
Return a pointer to the model.

Returns:
    ModelPtr: model setted.
        )", ( py::arg("self" )))
        .def("getMaterialField", &BehaviourProperty::getMaterialField, R"(
Return a pointer to the material field.

Returns:
    MaterialFieldPtr: material field setted.
        )", ( py::arg("self" )))
        .def("getBehaviourField", &BehaviourProperty::getBehaviourField, R"(
Return a pointer to the field for behaviour.

Returns:
    ConstantFieldOnCellsChar16Ptr: behaviour.
        )", ( py::arg("self" )))
        .def("getConvergenceCriteria", &BehaviourProperty::getConvergenceCriteria, R"(
Return a pointer to the field for convergence criteria.

Returns:
    ConstantFieldOnCellsRealPtr: convergence criteria.
        )", ( py::arg("self" )))
        .def("getMultipleBehaviourField", &BehaviourProperty::getMultipleBehaviourField, R"(
Return a pointer to the field for multiple behaviour like cristals.

Returns:
    ConstantFieldOnCellsChar16Ptr: multiple behaviour.
        )", ( py::arg("self" )));
};
