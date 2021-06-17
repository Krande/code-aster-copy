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

    bool ( MaterialProperty::*c1 )( std::string, const bool ) =
        &MaterialProperty::addPropertyReal;
    bool ( MaterialProperty::*c2 )( std::string, const ASTERDOUBLE &, const bool ) =
        &MaterialProperty::addPropertyReal;
    bool ( MaterialProperty::*c3 )( std::string, const bool ) =
        &MaterialProperty::addPropertyString;
    bool ( MaterialProperty::*c4 )( std::string, const std::string &, const bool ) =
        &MaterialProperty::addPropertyString;

    py::class_< MaterialProperty, MaterialPropertyPtr,
                py::bases< GenericMaterialProperty > >( "MaterialProperty", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< MaterialProperty, std::string >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< MaterialProperty, std::string, std::string >))
        .def( "addPropertyReal", c1 )
        .def( "addPropertyReal", c2 )
        .def( "addPropertyComplex", &MaterialProperty::addPropertyComplex )
        .def( "addPropertyString", c3 )
        .def( "addPropertyString", c4 )
        .def( "addPropertyFunction", &MaterialProperty::addPropertyFunction )
        .def( "addPropertyTable", &MaterialProperty::addPropertyTable )
        .def( "addPropertyVectorOfReal",
              &MaterialProperty::addPropertyVectorOfReal )
        .def( "addPropertyVectorOfFunction",
              &MaterialProperty::addPropertyVectorOfFunction )
        .def( "getName", &MaterialProperty::getName );
};
