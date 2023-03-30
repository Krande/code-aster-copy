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
        .def( "getNumberOfLines", &Table::getNumberOfLines )
        .def( "getParameters", &Table::getParameters )
        .def( "getColumnType", &Table::getColumnType )
        .def( "getColumn", &Table::getColumn );
    py::class_< TableOfFunctions, TableOfFunctions::TableOfFunctionsPtr, Table >(
        mod, "TableOfFunctions" )
        .def( py::init( &initFactoryPtr< TableOfFunctions > ) )
        .def( py::init( &initFactoryPtr< TableOfFunctions, std::string > ) )
        .def( "addFunction", &TableOfFunctions::addFunction )
        .def( "getFunction", &TableOfFunctions::getFunction )
        .def( "getNumberOfFunctions", &TableOfFunctions::getNumberOfFunctions );
};
