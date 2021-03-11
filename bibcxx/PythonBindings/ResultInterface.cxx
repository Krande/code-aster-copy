/**
 * @file ResultInterface.cxx
 * @brief Interface python de Result
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2020  EDF R&D                www.code-aster.org
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

#include "PythonBindings/ResultInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportResultToPython() {

    MaterialFieldPtr ( ResultClass::*c1 )() =
        &ResultClass::getMaterialField;
    MaterialFieldPtr ( ResultClass::*c2 )( int ) =
        &ResultClass::getMaterialField;

    ModelPtr ( ResultClass::*c3 )() =
        &ResultClass::getModel;
    ModelPtr ( ResultClass::*c4 )( int ) =
        &ResultClass::getModel;

    ElementaryCharacteristicsPtr ( ResultClass::*c5 )() =
        &ResultClass::getElementaryCharacteristics;
    ElementaryCharacteristicsPtr ( ResultClass::*c6 )( int ) =
        &ResultClass::getElementaryCharacteristics;

    bool ( ResultClass::*c7 )( const std::string ) const =
        &ResultClass::printMedFile;
    bool ( ResultClass::*c8 )( const std::string, std::string ) const =
        &ResultClass::printMedFile;

    py::class_< ResultClass, ResultClass::ResultPtr,
            py::bases< DataStructure > >( "Result", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< ResultClass, std::string >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< ResultClass, std::string, std::string >))
        .def( "addFieldOnNodesDescription", &ResultClass::addFieldOnNodesDescription )
        .def( "addMaterialField", &ResultClass::addMaterialField )
        .def( "addModel", &ResultClass::addModel )
        .def( "appendElementaryCharacteristicsOnAllRanks",
              &ResultClass::appendElementaryCharacteristicsOnAllRanks )
        .def( "appendMaterialFieldOnAllRanks",
              &ResultClass::appendMaterialFieldOnAllRanks )
        .def( "appendModelOnAllRanks", &ResultClass::appendModelOnAllRanks )
        .def( "listFields", &ResultClass::listFields )
        .def( "getElementaryCharacteristics", c5 )
        .def( "getElementaryCharacteristics", c6 )
        .def( "getMaterialField", c1 )
        .def( "getMaterialField", c2 )
        .def( "getMesh", &ResultClass::getMesh )
        .def( "getModel", c3 )
        .def( "getModel", c4 )
        .def( "getNumberOfRanks", &ResultClass::getNumberOfRanks )
        .def( "getAccessParameters", &ResultClass::getAccessParameters, R"(
Return the access parameters of the result as Python dict.

Returns:
    dict{str : list[int,float,str]}: Dict of values for each access variable.
        )", ( py::arg("self" )))
        .def( "getFieldsOnNodesNames", &ResultClass::getFieldsOnNodesNames, R"(
Return the names of the fields on nodes as Python list.

Returns:
    list[str]: List of names of the fields on nodes.
        )", ( py::arg("self" )))
        .def( "getFieldsOnCellsNames", &ResultClass::getFieldsOnCellsNames, R"(
Return the names of the fields on cells as Python list.

Returns:
    list[str]: List of names of the fields on cells.
        )", ( py::arg("self" )))
        .def( "getRanks", &ResultClass::getRanks )
        .def( "getRealFieldOnNodes", &ResultClass::getRealFieldOnNodes )
        .def( "getRealFieldOnCells", &ResultClass::getRealFieldOnCells )
        .def( "printMedFile", c7 )
        .def( "printMedFile", c8 )
        .def( "setMesh", &ResultClass::setMesh )
        .def( "update", &ResultClass::update )

        .def( "getTable", &ListOfTablesClass::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) )
        ;
};
