/**
 * @file ElementaryMatrixInterface.cxx
 * @brief Interface python de ElementaryMatrix
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

/* person_in_charge: nicolas.sellenet at edf.fr */

#include <boost/python.hpp>

namespace py = boost::python;
#include "PythonBindings/ElementaryTermInterface.h"

#include <PythonBindings/factory.h>

void exportElementaryTermToPython() {

    py::class_< ElementaryTermReal, ElementaryTermRealPtr, py::bases< DataField > >(
        "ElementaryTermReal", py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< ElementaryTermReal > ) )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< ElementaryTermReal, std::string > ) )
        .def( "getFiniteElementDescriptor", &ElementaryTermReal::getFiniteElementDescriptor, R"(
Return the finite element descriptor

Returns:
    FiniteElementDescriptor: finite element descriptor
        )",
              ( py::arg( "self" ) ) )
        .def( "getOption", &ElementaryTermReal::getOption, R"(
Return the optior used to compute it

Returns:
    str: name of the option
        )",
              ( py::arg( "self" ) ) )
        .def( "getMesh", &ElementaryTermReal::getMesh, R"(
Return the mesh

Returns:
    BaseMesh: a pointer to the mesh
        )",
              ( py::arg( "self" ) ) )
        .def( "getPhysicalQuantity", &ElementaryTermReal::getPhysicalQuantity, R"(
Return the physical quantity

Returns:
    str: name of the physical quantity
        )",
              ( py::arg( "self" ) ) );

    py::class_< ElementaryTermComplex, ElementaryTermComplexPtr, py::bases< DataField > >(
        "ElementaryTermComplex", py::no_init )
        .def( "__init__", py::make_constructor( &initFactoryPtr< ElementaryTermComplex > ) )
        .def( "__init__",
              py::make_constructor( &initFactoryPtr< ElementaryTermComplex, std::string > ) )
        .def( "getFiniteElementDescriptor", &ElementaryTermComplex::getFiniteElementDescriptor, R"(
Return the finite element descriptor

Returns:
    FiniteElementDescriptor: finite element descriptor
        )",
              ( py::arg( "self" ) ) )
        .def( "getOption", &ElementaryTermComplex::getOption, R"(
Return the optior used to compute it

Returns:
    str: name of the option
        )",
              ( py::arg( "self" ) ) )
        .def( "getMesh", &ElementaryTermComplex::getMesh, R"(
Return the mesh

Returns:
    BaseMesh: a pointer to the mesh
        )",
              ( py::arg( "self" ) ) )
        .def( "getPhysicalQuantity", &ElementaryTermComplex::getPhysicalQuantity, R"(
Return the physical quantity

Returns:
    str: name of the physical quantity
        )",
              ( py::arg( "self" ) ) );
};
