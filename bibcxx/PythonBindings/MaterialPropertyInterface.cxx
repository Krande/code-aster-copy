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


#include "PythonBindings/MaterialPropertyInterface.h"
#include "PythonBindings/BaseMaterialPropertyInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportMaterialPropertyToPython() {

    bool ( MaterialPropertyClass::*c1 )( std::string, const bool ) =
        &MaterialPropertyClass::addRealProperty;
    bool ( MaterialPropertyClass::*c2 )( std::string, const ASTERDOUBLE &, const bool ) =
        &MaterialPropertyClass::addRealProperty;
    bool ( MaterialPropertyClass::*c3 )( std::string, const bool ) =
        &MaterialPropertyClass::addStringProperty;
    bool ( MaterialPropertyClass::*c4 )( std::string, const std::string &, const bool ) =
        &MaterialPropertyClass::addStringProperty;

    py::class_< MaterialPropertyClass, MaterialPropertyPtr,
                py::bases< GenericMaterialPropertyClass > >( "MaterialProperty", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< MaterialPropertyClass, std::string >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< MaterialPropertyClass, std::string, std::string >))
        .def( "addRealProperty", c1 )
        .def( "addRealProperty", c2 )
        .def( "addComplexProperty", &MaterialPropertyClass::addComplexProperty )
        .def( "addStringProperty", c3 )
        .def( "addStringProperty", c4 )
        .def( "addFunctionProperty", &MaterialPropertyClass::addFunctionProperty )
        .def( "addTableProperty", &MaterialPropertyClass::addTableProperty )
        .def( "addVectorOfRealProperty",
              &MaterialPropertyClass::addVectorOfRealProperty )
        .def( "addVectorOfFunctionProperty",
              &MaterialPropertyClass::addVectorOfFunctionProperty )
        .def( "getName", &MaterialPropertyClass::getName );
};
