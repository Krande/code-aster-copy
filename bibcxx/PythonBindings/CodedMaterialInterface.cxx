/**
 * @file CodedMaterialInterface.cxx
 * @brief Interface python de CodedMaterial
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

#include "PythonBindings/CodedMaterialInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS( exec_allocate, allocate, 0, 1 )

void exportCodedMaterialToPython() {

    py::class_< CodedMaterial, CodedMaterialPtr > c1(
        "CodedMaterial", py::no_init );
    c1.def( "__init__",
            py::make_constructor(
                &initFactoryPtr< CodedMaterial, MaterialFieldPtr, ModelPtr >));
    c1.def( "__init__",
            py::make_constructor(&initFactoryPtr< CodedMaterial, std::string,
                                              MaterialFieldPtr, ModelPtr >));
    c1.def( "allocate", &CodedMaterial::allocate, exec_allocate() );
    c1.def( "constant", &CodedMaterial::constant );
};
