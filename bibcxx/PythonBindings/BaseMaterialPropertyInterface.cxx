/**
 * @file MaterialPropertyInterface.cxx
 * @brief Interface python de MaterialProperty
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

#include "PythonBindings/BaseMaterialPropertyInterface.h"

#include "aster_pybind.h"

void exportBaseMaterialPropertyToPython( py::module_ &mod ) {

    py::class_< GenericMaterialProperty, GenericMaterialPropertyPtr >( mod,
                                                                       "GenericMaterialProperty" )
        // fake initFactoryPtr: created by subclasses
        .def( py::init( &initFactoryPtr< GenericMaterialProperty > ) )
        .def( "getAsterName", &GenericMaterialProperty::getAsterName )
        .def_property( "hasTractionFunction", &GenericMaterialProperty::hasTractionFunction,
                       &GenericMaterialProperty::setHasTractionFunction, R"(
bool: Attribute that holds the need of a traction function.
                       )" )
        .def( "getValueComplex", &GenericMaterialProperty::getValueComplex )
        .def( "getValueReal", &GenericMaterialProperty::getValueReal )
        .def( "getValueString", &GenericMaterialProperty::getValueString )
        .def( "getValueGenericFunction", &GenericMaterialProperty::getValueGenericFunction )
        .def( "getNumberOfListOfPropertiesReal",
              &GenericMaterialProperty::getNumberOfListOfPropertiesReal )
        .def( "getNumberOfListOfPropertiesFunction",
              &GenericMaterialProperty::getNumberOfListOfPropertiesFunction )
        .def( "getValueTable", &GenericMaterialProperty::getValueTable )
        .def( "hasValueComplex", &GenericMaterialProperty::hasValueComplex )
        .def( "hasValueReal", &GenericMaterialProperty::hasValueReal )
        .def( "hasValueString", &GenericMaterialProperty::hasValueString )
        .def( "hasValueGenericFunction", &GenericMaterialProperty::hasValueGenericFunction )
        .def( "hasValueTable", &GenericMaterialProperty::hasValueTable )
        // functions will be tried in sequence, order is important to avoid unexpected implicit
        // conversions
        .def( "setValue",
              py::overload_cast< std::string, ASTERDOUBLE >( &GenericMaterialProperty::setValue ) )
        .def( "setValue",
              py::overload_cast< std::string, ASTERCOMPLEX >( &GenericMaterialProperty::setValue ) )
        .def( "setValue",
              py::overload_cast< std::string, std::string >( &GenericMaterialProperty::setValue ) )
        .def( "setValue",
              py::overload_cast< std::string, FunctionPtr >( &GenericMaterialProperty::setValue ) )
        .def( "setValue",
              py::overload_cast< std::string, TablePtr >( &GenericMaterialProperty::setValue ) )
        .def( "setValue", py::overload_cast< std::string, Function2DPtr >(
                              &GenericMaterialProperty::setValue ) )
        .def( "setValue",
              py::overload_cast< std::string, VectorReal >( &GenericMaterialProperty::setValue ) )
        .def( "setValue", py::overload_cast< std::string, VectorFunction >(
                              &GenericMaterialProperty::setValue ) )
        .def( "setValue",
              py::overload_cast< std::string, FormulaPtr >( &GenericMaterialProperty::setValue ) )

        .def( "setSortedListParameters", &GenericMaterialProperty::setSortedListParameters );
};
