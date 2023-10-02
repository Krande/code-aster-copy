/**
 * @file TableInterface.cxx
 * @brief Interface python de Table
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
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

#include "PythonBindings/TableInterface.h"

#include "aster_pybind.h"

#include "PythonBindings/DataStructureInterface.h"

void exportTableToPython( py::module_ &mod ) {

    py::class_< Table, Table::TablePtr, DataStructure >( mod, "Table" )
        .def( py::init( &initFactoryPtr< Table > ) )
        .def( py::init( &initFactoryPtr< Table, std::string > ) )
        .def( "getNumberOfLines", &Table::getNumberOfLines, R"(
Returns the number of lines of the table.

Returns:
    int: Number of lines.
        )" )
        .def( "getParameters", &Table::getParameters, R"(
Return the parameters names.

Returns:
    list[str]: Names of the parameters.
        )" )
        .def( "getColumnType", &Table::getColumnType, R"(
Return the type of values in a column.

Arguments:
    param (str): Parameter name.

Returns:
    str: "I" for integers, "R" for reals, "C" for complex, "Knn" for strings.
        )",
              py::arg( "param" ) )
        .def( "getValues", &Table::getValues, R"(
For internal use only. See *get_column()*.
        )" );
    py::class_< TableOfFunctions, TableOfFunctions::TableOfFunctionsPtr, Table >(
        mod, "TableOfFunctions" )
        .def( py::init( &initFactoryPtr< TableOfFunctions > ) )
        .def( py::init( &initFactoryPtr< TableOfFunctions, std::string > ) )
        .def( "addFunction", &TableOfFunctions::addFunction, R"(
Add a function into the table.
        )" )
        .def( "getFunction", &TableOfFunctions::getFunction, R"(
Returns the function stored at a given position.

Arguments:
    pos [int]: Index of the function to return (0-based).

Returns:
    *Function*: Function stored.
        )",
              py::arg( "pos" ) )
        .def( "getNumberOfFunctions", &TableOfFunctions::getNumberOfFunctions, R"(
Returns the number of functions stored in the table.

Returns:
    int: Number of functions.
        )" );
};
