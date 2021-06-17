/**
 * @file MechanicalLoadInterface.cxx
 * @brief Interface python de MechanicalLoad
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
#include "PythonBindings/MechanicalLoadInterface.h"

void exportMechanicalLoadToPython() {


    py::class_< MechanicalLoadReal,
                MechanicalLoadReal::MechanicalLoadPtr,
                py::bases< DataStructure > >( "MechanicalLoadReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< MechanicalLoadReal, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< MechanicalLoadReal, std::string, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              &MechanicalLoadReal::getFiniteElementDescriptor )
        .def( "hasLoad", &MechanicalLoadReal::hasLoad )
        .def( "updateValuePointers", &MechanicalLoadReal::updateValuePointers)
        .def( "getModel", &MechanicalLoadReal::getModel)
        .def( "getMesh", &MechanicalLoadReal::getMesh)
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) );

    py::class_<  MechanicalLoadFunction,
                 MechanicalLoadFunction::MechanicalLoadPtr,
                py::bases< DataStructure > >( "MechanicalLoadFunction", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr<  MechanicalLoadFunction, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr<  MechanicalLoadFunction, std::string, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              & MechanicalLoadFunction::getFiniteElementDescriptor )
        .def( "hasLoad", &MechanicalLoadFunction::hasLoad )
        .def( "updateValuePointers", &MechanicalLoadFunction::updateValuePointers)
        .def( "getModel", & MechanicalLoadFunction::getModel)
        .def( "getMesh", & MechanicalLoadFunction::getMesh)
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) );


    py::class_<  MechanicalLoadComplex,
                 MechanicalLoadComplex::MechanicalLoadPtr,
                py::bases< DataStructure > >( "MechanicalLoadComplex", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr<  MechanicalLoadComplex, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr<  MechanicalLoadComplex, std::string, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              & MechanicalLoadComplex::getFiniteElementDescriptor )
        .def( "hasLoad", &MechanicalLoadComplex::hasLoad )
        .def( "updateValuePointers", &MechanicalLoadComplex::updateValuePointers)
        .def( "getModel", & MechanicalLoadComplex::getModel)
        .def( "getMesh", & MechanicalLoadComplex::getMesh)
        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) );

};
