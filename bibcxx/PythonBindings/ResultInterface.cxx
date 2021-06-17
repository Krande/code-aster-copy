/**
 * @file ResultInterface.cxx
 * @brief Interface python de Result
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

#include "PythonBindings/ResultInterface.h"
#include "PythonBindings/factory.h"
#include <boost/python.hpp>

namespace py = boost::python;

void exportResultToPython() {

    MaterialFieldPtr ( Result::*c1 )() =
        &Result::getMaterialField;
    MaterialFieldPtr ( Result::*c2 )( int ) =
        &Result::getMaterialField;

    ModelPtr ( Result::*c3 )() =
        &Result::getModel;
    ModelPtr ( Result::*c4 )( int ) =
        &Result::getModel;

    ElementaryCharacteristicsPtr ( Result::*c5 )() =
        &Result::getElementaryCharacteristics;
    ElementaryCharacteristicsPtr ( Result::*c6 )( int ) =
        &Result::getElementaryCharacteristics;

    bool ( Result::*c7 )( const std::string ) const =
        &Result::printMedFile;
    bool ( Result::*c8 )( const std::string, std::string ) const =
        &Result::printMedFile;

    py::class_< Result, Result::ResultPtr,
            py::bases< DataStructure > >( "Result", py::no_init )
        .def( "__init__",
              py::make_constructor(&initFactoryPtr< Result, std::string >))
        .def( "__init__",
              py::make_constructor(
                  &initFactoryPtr< Result, std::string, std::string >))
        .def( "addFieldOnNodesDescription", &Result::addFieldOnNodesDescription )
        .def( "addMaterialField", &Result::addMaterialField )
        .def( "addModel", &Result::addModel )
        .def( "addElementaryCharacteristics", &Result::addElementaryCharacteristics )
        .def( "appendElementaryCharacteristicsOnAllRanks",
              &Result::appendElementaryCharacteristicsOnAllRanks )
        .def( "appendMaterialFieldOnAllRanks",
              &Result::appendMaterialFieldOnAllRanks )
        .def( "appendModelOnAllRanks", &Result::appendModelOnAllRanks )
        .def( "listFields", &Result::listFields )
        .def( "getAllElementaryCharacteristics", &Result::getAllElementaryCharacteristics, R"(
Return the list of all elementary characteristics used in the result

Returns:
    list[ElementaryCharacteristicsPtr]: list of ElementaryCharacteristics.
        )", ( py::arg("self" )) )
        .def( "getElementaryCharacteristics", c5 )
        .def( "getElementaryCharacteristics", c6 )
        .def( "getMaterialFields", &Result::getMaterialFields, R"(
Return the list of all material fields used in the result

Returns:
    list[MaterialFieldPtr]: list of material field.
        )", ( py::arg("self" )) )
        .def( "getMaterialField", c1 )
        .def( "getMaterialField", c2 )
        .def( "getMesh", &Result::getMesh )
        .def( "getModels", &Result::getModels, R"(
Return the list of all models used in the result

Returns:
    list[ModelPtr]: list of models.
        )", ( py::arg("self" )) )
        .def( "getModel", c3 )
        .def( "getModel", c4 )
        .def( "getNumberOfRanks", &Result::getNumberOfRanks )
        .def( "getAccessParameters", &Result::getAccessParameters, R"(
Return the access parameters of the result as Python dict.

Returns:
    dict{str : list[int,float,str]}: Dict of values for each access variable.
        )", ( py::arg("self" )))
        .def( "getFieldsOnNodesNames", &Result::getFieldsOnNodesNames, R"(
Return the names of the fields on nodes as Python list.

Returns:
    list[str]: List of names of the fields on nodes.
        )", ( py::arg("self" )))
        .def( "getFieldsOnCellsNames", &Result::getFieldsOnCellsNames, R"(
Return the names of the fields on cells as Python list.

Returns:
    list[str]: List of names of the fields on cells.
        )", ( py::arg("self" )))
        .def( "getRanks", &Result::getRanks )
        .def( "getFieldOnNodesReal", &Result::getFieldOnNodesReal )
        .def( "getFieldOnCellsReal", &Result::getFieldOnCellsReal )
        .def( "printMedFile", c7 )
        .def( "printMedFile", c8 )
        .def( "setMesh", &Result::setMesh )
        .def( "build", &Result::build )
        .def( "printInfo", &Result::printInfo )

        .def( "getTable", &ListOfTables::getTable, R"(
Extract a Table from the datastructure.

Arguments:
    identifier (str): Table identifier.

Returns:
    Table: Table stored with the given identifier.
        )",
              ( py::arg( "self" ), py::arg( "identifier" ) ) )
        ;
};
