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


    py::class_< MechanicalLoadRealClass,
                MechanicalLoadRealClass::MechanicalLoadPtr,
                py::bases< DataStructure > >( "MechanicalLoadReal", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< MechanicalLoadRealClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< MechanicalLoadRealClass, std::string, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              &MechanicalLoadRealClass::getFiniteElementDescriptor )
        .def( "hasLoad", &MechanicalLoadRealClass::hasLoad )
        .def( "updateValuePointers", &MechanicalLoadRealClass::updateValuePointers)
        .def( "getModel", &MechanicalLoadRealClass::getModel)
        .def( "getMesh", &MechanicalLoadRealClass::getMesh)
        .def( "getTable", &ListOfTablesClass::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) );

    py::class_<  MechanicalLoadFunctionClass,
                 MechanicalLoadFunctionClass::MechanicalLoadPtr,
                py::bases< DataStructure > >( "MechanicalLoadFunction", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr<  MechanicalLoadFunctionClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr<  MechanicalLoadFunctionClass, std::string, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              & MechanicalLoadFunctionClass::getFiniteElementDescriptor )
        .def( "hasLoad", &MechanicalLoadFunctionClass::hasLoad )
        .def( "updateValuePointers", &MechanicalLoadFunctionClass::updateValuePointers)
        .def( "getModel", & MechanicalLoadFunctionClass::getModel)
        .def( "getMesh", & MechanicalLoadFunctionClass::getMesh)
        .def( "getTable", &ListOfTablesClass::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) );


    py::class_<  MechanicalLoadComplexClass,
                 MechanicalLoadComplexClass::MechanicalLoadPtr,
                py::bases< DataStructure > >( "MechanicalLoadComplex", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr<  MechanicalLoadComplexClass, ModelPtr >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr<  MechanicalLoadComplexClass, std::string, ModelPtr >))
        .def( "getFiniteElementDescriptor",
              & MechanicalLoadComplexClass::getFiniteElementDescriptor )
        .def( "hasLoad", &MechanicalLoadComplexClass::hasLoad )
        .def( "updateValuePointers", &MechanicalLoadComplexClass::updateValuePointers)
        .def( "getModel", & MechanicalLoadComplexClass::getModel)
        .def( "getMesh", & MechanicalLoadComplexClass::getMesh)
        .def( "getTable", &ListOfTablesClass::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) );

};
