/**
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

#include "PythonBindings/MGISBehaviourInterface.h"

#include "aster_pybind.h"

void exportMGISBehaviourToPython( py::module_ &mod ) {

#ifdef ASTER_HAVE_MGIS
    py::class_< MGISBehaviour, MGISBehaviour::MGISBehaviourPtr, DataStructure >( mod,
                                                                                 "MGISBehaviour" )
        .def( py::init( &initFactoryPtr< MGISBehaviour > ) )
        .def( py::init( &initFactoryPtr< MGISBehaviour, std::string > ) )
        .def( "setLibPath", &MGISBehaviour::setLibPath, R"(
Set the path to the MFront library.

Arguments:
    path: Library path.
        )",
              py::arg( "path" ) )
        .def( "setBehaviourName", &MGISBehaviour::setBehaviourName, R"(
Define the name of the behaviour to be used from the MFront library.

Arguments:
    name: Name of the behaviour.
        )",
              py::arg( "name" ) )

        .def( "_internal_state", &MGISBehaviour::internal_state, R"(
Returns internal state.

Returns:
    list[str]: internal state.
        )" );
#endif
};
