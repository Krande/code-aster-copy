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

#include "PythonBindings/MaterialPropertyInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/BaseMaterialPropertyInterface.h"

void exportMaterialPropertyToPython( py::module_ &mod ) {

    py::class_< MaterialProperty, MaterialPropertyPtr, GenericMaterialProperty >(
        mod, "MaterialProperty" )
        .def( py::init( &initFactoryPtr< MaterialProperty, std::string > ) )
        .def( py::init( &initFactoryPtr< MaterialProperty, std::string, std::string > ) )
        .def( "addPropertyReal",
              py::overload_cast< std::string, const bool >( &MaterialProperty::addPropertyReal ) )
        .def( "addPropertyReal", py::overload_cast< std::string, const ASTERDOUBLE &, const bool >(
                                     &MaterialProperty::addPropertyReal ) )
        .def( "addPropertyComplex", &MaterialProperty::addPropertyComplex )
        .def( "addPropertyString",
              py::overload_cast< std::string, const bool >( &MaterialProperty::addPropertyString ) )
        .def( "addPropertyString",
              py::overload_cast< std::string, const std::string &, const bool >(
                  &MaterialProperty::addPropertyString ) )
        .def( "addPropertyFunction", &MaterialProperty::addPropertyFunction )
        .def( "addPropertyTable", &MaterialProperty::addPropertyTable )
        .def( "addPropertyVectorOfReal", &MaterialProperty::addPropertyVectorOfReal )
        .def( "addPropertyVectorOfFunction", &MaterialProperty::addPropertyVectorOfFunction )
        .def( "getName", &MaterialProperty::getName );
};
