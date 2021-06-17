/**
 * @file MaterialPropertyInterface.cxx
 * @brief Interface python de MaterialProperty
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

#include "PythonBindings/BaseMaterialPropertyInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportBaseMaterialPropertyToPython() {

    bool ( GenericMaterialProperty::*c1 )(  ) const =
        &GenericMaterialProperty::hasTractionFunction;
    void ( GenericMaterialProperty::*c2 )( const bool )=
        &GenericMaterialProperty::hasTractionFunction;

    // Candidates for setValue
    bool ( GenericMaterialProperty::*d1 )( std::string, ASTERDOUBLE ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d2 )( std::string, ASTERCOMPLEX ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d3 )( std::string, std::string ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d4 )( std::string, FunctionPtr ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d5 )( std::string, TablePtr ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d6 )( std::string, Function2DPtr ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d7 )( std::string, VectorReal ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d8 )( std::string, VectorFunction ) =
        &GenericMaterialProperty::setValue;
    bool ( GenericMaterialProperty::*d9 )( std::string, FormulaPtr ) =
        &GenericMaterialProperty::setValue;

    py::class_< GenericMaterialProperty, GenericMaterialPropertyPtr >(
        "GenericMaterialProperty", py::no_init )
        // fake initFactoryPtr: created by subclasses
        .def( "__init__", py::make_constructor(&initFactoryPtr< GenericMaterialProperty >))
        .def( "getAsterName", &GenericMaterialProperty::getAsterName )
        .def( "hasTractionFunction", c1 )
        .def( "hasTractionFunction", c2 )
        .def( "getValueComplex", &GenericMaterialProperty::getValueComplex )
        .def( "getValueReal", &GenericMaterialProperty::getValueReal )
        .def( "getValueString", &GenericMaterialProperty::getValueString )
        .def( "getValueGenericFunction",
              &GenericMaterialProperty::getValueGenericFunction )
        .def( "getNumberOfListOfPropertiesReal",
              &GenericMaterialProperty::getNumberOfListOfPropertiesReal )
        .def( "getNumberOfListOfPropertiesFunction",
              &GenericMaterialProperty::getNumberOfListOfPropertiesFunction )
        .def( "getValueTable", &GenericMaterialProperty::getValueTable )
        .def( "hasValueComplex", &GenericMaterialProperty::hasValueComplex )
        .def( "hasValueReal", &GenericMaterialProperty::hasValueReal )
        .def( "hasValueString", &GenericMaterialProperty::hasValueString )
        .def( "hasValueGenericFunction",
              &GenericMaterialProperty::hasValueGenericFunction )
        .def( "hasValueTable", &GenericMaterialProperty::hasValueTable )
        // J'ai changé le nom pour setValueReal car il y a une conversion implite
        // float -> complex sinon (à changer un jour)
        .def( "setValueReal", d1 )
        .def( "setValue", d2 )
        .def( "setValue", d3 )
        .def( "setValue", d4 )
        .def( "setValue", d5 )
        .def( "setValue", d6 )
        .def( "setValueVectorReal", d7 )
        .def( "setValue", d8 )
        .def( "setValue", d9 )
        .def( "setSortedListParameters",
              &GenericMaterialProperty::setSortedListParameters );

};
