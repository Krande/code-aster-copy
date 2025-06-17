/**
 * @file DataStructureInterface.cxx
 * @brief Interface python de DataStructure
 * @author Nicolas Sellenet
 * @section LICENCE
 *   Copyright (C) 1991 - 2025  EDF R&D                www.code-aster.org
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

#include "PythonBindings/DataStructureInterface.h"

#include "aster_pybind.h"

#include <iostream>

void exportDataStructureToPython( py::module_ &mod ) {

    py::class_< DataStructure, DataStructurePtr >( mod, "DataStructure" )
        // fake initFactoryPtr: created by subclasses
        // fake initFactoryPtr: created by subclasses
        .def_property( "_ptr_sdj", &DataStructure::getSDJ, &DataStructure::setSDJ )
        .def_property( "_ptr_cache", &DataStructure::getCache, &DataStructure::setCache )
        .def( "id", &DataStructure::id,
              R"(
Return the identity of the object.

Returns:
    int: Identifier (address as int).
        )" )

        .def( "addDependency", &DataStructure::addDependency,
              R"(
Add a dependency to a *DataStructure*.

Arguments:
    ds (*DataStructure*): Parent *DataStructure* to depend on.
        )",
              py::arg( "ds" ) )

        .def( "removeDependency", &DataStructure::removeDependency,
              R"(
Remove a dependency to a *DataStructure*.

Arguments:
    ds (*DataStructure*): Parent *DataStructure* to be removed from
        dependencies.
        )",
              py::arg( "ds" ) )

        .def( "resetDependencies", &DataStructure::resetDependencies,
              R"(
Clear the list of explicit dependencies.
        )" )

        .def( "getDependencies", &DataStructure::getDependencies,
              R"(
Return the explicit dependencies.

Returns:
    list[*DataStructure*]: List of parents (dependencies) *DataStructure*.
        )" )

        .def( "getName", &DataStructure::getName,
              R"(
Return the internal (*Jeveux*) name of the *DataStructure*.

Returns:
    str: Internal/*Jeveux* name.
        )" )

        .def_property( "userName", &DataStructure::getUserName, &DataStructure::setUserName,
                       py::return_value_policy::copy,
                       R"(
str: Name of the user variable that holds this object.
        )" )
        .def( "getType", &DataStructure::getType,
              R"(
Return the name of the *DataStructure* type.

Returns:
    str: Name of the *DataStructure* type.
        )" )
        .def( "setTitle", &DataStructure::setTitle,
              R"(
Set the tile of the *DataStructure* .

Arguments:
    title [str]: Title of the *DataStructure*.
        )",
              py::arg( "title" ) )
        .def( "getTitle", &DataStructure::getTitle,
              R"(
Return the tile of the *DataStructure* .

Returns:
    str: Title of the *DataStructure*.
        )" )
        .def( "debugPrint",
              py::overload_cast< int, bool >( &DataStructure::debugPrint, py::const_ ), R"(
Print the raw content of a *DataStructure* on the selected file.

Args:
    unit (int): File number (default: 6, means stdout).
    synchro (bool): To synchronize prints between processors (default: True).
        )",
              py::arg( "unit" ) = 6, py::arg( "synchro" ) = true )
        .def( "build", &DataStructure::build,
              R"(
Update the *DataStructure* attributes from the *Jeveux* objects.
*Only use internally after calling fortran subroutines*.

Returns:
    bool: *True* if all went ok, *False* otherwise.
        )" );

    py::class_< DSWithCppPickling, DSWithCppPicklingPtr, DataStructure >( mod,
                                                                          "DSWithCppPickling" );
    // fake initFactoryPtr: created by subclasses
    // fake initFactoryPtr: created by subclasses
};
