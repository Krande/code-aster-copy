/**
 * @file ListOfLoadsInterface.cxx
 * @brief Interface python de ListOfLoads
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

#include "PythonBindings/ListOfLoadsInterface.h"

#include "PythonBindings/LoadUtilities.h"
#include "PythonBindings/factory.h"

#include <boost/python.hpp>

namespace py = boost::python;

void exportListOfLoadsToPython() {

    py::class_< ListOfLoads, ListOfLoadsPtr, py::bases< DataStructure > > c1( "ListOfLoads",
                                                                              py::no_init );
    c1.def( "__init__", py::make_constructor( &initFactoryPtr< ListOfLoads > ) );
    c1.def( "__init__", py::make_constructor( &initFactoryPtr< ListOfLoads, std::string > ) );
    c1.def( "__init__", py::make_constructor( &initFactoryPtr< ListOfLoads, ModelPtr > ) );
    c1.def( "__init__",
            py::make_constructor( &initFactoryPtr< ListOfLoads, std::string, ModelPtr > ) );

    c1.def( "getModel", &ListOfLoads::getModel, R"(
Return the model used

Returns:
    Model: model used
        )",
            ( py::arg( "self" ) ) );
};
